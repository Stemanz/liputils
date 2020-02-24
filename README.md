# liputils
A small Python package to manipulate complex lipids.

## Overview
```liputils``` makes it easy to strip fatty acids-like residues from individual molecular lipids. This is done by ```liputils``` by reading the lipid name.

Tracking individual residues is is particularly useful when wanting to track how the carbon chains move across the lipidome, independently from where they are attached to. For instance, it is possible to see if the general trend of long carbon residues in the plasma matches the data available from the dietary treatment.

The ```Lipid``` class takes care of extracting information from the lipid name:

```python
>>> from liputils import Lipid

>>> l = Lipid("PG 18:1/20:1", amount=0.012512)

>>> l.mass
802.5724

>>> l.lipid_class()
'PG'

>>> l.name
'PG 18:1/20:1'

>>> l.residues()
(['18:1', '20:1'], 1)

>>> l.molecules
7534902640.2784
```
\
The number of molecules is calculated from the ```amount``` parameter, defaulting to picomoles. This can be changed:

```python
>>> l = Lipid("PG 18:1/20:1", amount=0.012512, unit="femtomoles")

>>> l.molecules
7534902.640278401
```
\
In the case of unresolved ambiguities of the lipid isomers, it is possible to either extract all of them and choose how to manage that information by taking into consideration how many ambiguities there are:

```python
>>> l = Lipid("TAG 48:2 total (14:0/16:0/18:2)(14:0/16:1/18:1)(16:0/16:1/16:1)")

>>> l.residues()                
(['14:0', '16:0', '18:2', '14:0', '16:1', '18:1', '16:0', '16:1', '16:1'], 3)
```
\
Or, it is possible to reject non unambiguous lipids altogether by calling ```.residues()``` with the ```drop_ambiguous``` parameter:

```python
>>> l = Lipid("PG 18:1/20:1")               

>>> l.residues(drop_ambiguous=True)             
(['18:1', '20:1'], 1)


>>> l = Lipid("TAG 48:2 total (14:0/16:0/18:2)(14:0/16:1/18:1)(16:0/16:1/16:1)")

>>> l.residues(drop_ambiguous=True)          
([], 0)
```

## One-step lipidomics data conversion
Lipidomics data should be loaded in a ```pandas.DataFrame``` table. The accepted format is a vertical index with lipid names, and samples in column. Just like this:

![](https://github.com/Stemanz/liputils/raw/master/images/liputils_sample_table.png)


```make_residues_table``` will take care of dropping non-numerical columns, as well as to trim the lipid list of elements that should not be processed, like total lipid class counts. These can be further specified through the ```unwanted``` parameter.
Getting the transformed table is super easy:
```python
# df is out dataframe
res = make_residues_table(df)
```
\
In ```res```, we will find the resulting table:

![](https://github.com/Stemanz/liputils/raw/master/images/liputils_processed_sample_table.png)


That's it! For further info, don't forget to investigate around:

```python
help(make_residues_table)
```

```
Parameters
==========

dataframe: a pandas dataframe of data. Lipid names as index, and samples as columns
    (just unlike sklearn wants it, but as you might get it from Tableau software
    tables. Just dataframe.T your table - that would just do the trick).

drop_ambiguous: <bool> don't take isobars into consideration. Defaults to False. If True,
    each residue is divided by its uncertainty.

name: <str> a tag that gets attached to the returned dataframe, so you can use it
    to save it afterwards. The tag is found in the .name attribute.

replace_nan: <object> the object you would like to replace your missing values with.
    It can be set to False, but I would suggest against what.

cleanup: <bool> Whether to perform a cleanup of unwanted lipids that can be present
    in the index. Unwanted strings are read from the 'unwanted' parameter. Defaults
    to True

absolute_amount <bool> Wheter to count the individual number of residues, rather to
    sticking to the same units found in the original table. Defaults to False

unwanted: <list> <set> <tuple> Strings that must be removed from the lipid index. Defaults
    to ["total", "fc", "tc"]

returns:
========

pandas DataFrame
```
