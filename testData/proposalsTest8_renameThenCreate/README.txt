Creating a taxon in an already moved/renamed taxa should be a warning, not an error. 

Here, the last line of the sheet renames the class Ellioviricetes to Bunyaviricetes, then the other lines, which are above, create things, like the order Goldbavirales, in Ellioviricetes. 

However, the processor does higher-ranked changes, first, so while it recognized the fact that it had already been renamed, it threw an error, instead of going with the flow. 

This sheet will still toss a lot of CREATE.PARENT_LINEAGE warnings, which they could avoid by creating the new taxa in the new order, 'Ellioviricetes'.
