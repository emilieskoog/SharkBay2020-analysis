# README

#### This will walk you through a simple tutorial on how I used the folium package in python to create a map of Shark Bay, Australia.

>Note: For more information on how to use the folium package, visit [this page](https://python-visualization.github.io/folium/). 


## :memo: Mapping Shark Bay

### Step 1: Import packages

Open python and import packages.

```
import folium
import pandas as pd
```

### Step 2: Create and view map

I set my location as -26.141967, 114.236685. 

```
Shark_Bay_map = folium.Map([-26.141967, 114.236685], zoom_start=8, width=600, height=500, tiles='Stamen Watercolor', control_scale=True)
Shark_Bay_map
```
Executing this command should give you something like this:


![](https://i.imgur.com/ukjimEF.jpg)

### Step 3: Save your map!

```
Shark_Bay_map.save('/insert/your/desired/path/Shark_Bay_map.html')
```
