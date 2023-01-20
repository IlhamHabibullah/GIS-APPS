#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ee

print(ee.__version__)


# In[2]:


import ee
ee.Authenticate()


# In[4]:


ee.Initialize()


# In[5]:


Map = geemap.Map()
Map


# In[8]:


aoi = ee.Geometry.Polygon(
        [[[109.10315235495382, -0.49113846815019824],
          [109.10315235495382, -0.6332672093911977],
          [109.3324919545632, -0.6332672093911977],
          [109.3324919545632, -0.49113846815019824]]], None, False)

#
# Script modified by Ilham Habibullah - Kelautan
#

#################################/
############/   COMPUTE MVI   ############/
#################################/

#=================================================================
# STEP 1. CREATE FILTERS

# Center the map on AOI
bound = aoi.bounds()
Map.centerObject(bound, 12)

#Construct start and end dates:
start = ee.Date('2020-01-01')
finish = ee.Date('2020-12-31')

#=================================================================
# STEP 2. LOAD LANDSAT 8 IMAGE COLLECTION

# Load Landsat 8 surface reflectance data
l8sr = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')             .filterBounds(aoi)             .filterDate(start, finish)

#=================================================================
# STEP 3. CREATE A CLOUD FREE MOSAIC

def maskL8sr(image):
  cloudShadowBitMask = ee.Number(2).pow(3).int()
  cloudsBitMask = ee.Number(2).pow(5).int()
  qa = image.select('pixel_qa')
  mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)       .And(qa.bitwiseAnd(cloudsBitMask).eq(0))
  return image.updateMask(mask).divide(10000)

composite = l8sr.map(maskL8sr)                     .reduce(ee.Reducer.median())

#=================================================================
# STEP 4. MASKING BEFORE ANALYSIS

# Masking for pixel above 50 m
srtm = ee.Image('USGS/SRTMGL1_003')
elevation = srtm.select('elevation')
masksrtm = srtm.lt(50)
maskedsrtm = composite.updateMask(masksrtm)

# Water masking
hansenImage = ee.Image('UMD/hansen/global_forest_change_2015')
datamask = hansenImage.select('datamask')
maskland = datamask.eq(1)
maskedcomposite = maskedsrtm.updateMask(maskland)
lcomposite = maskedcomposite.clip(aoi)

#=================================================================
# STEP 5. COMPUTE THE MVI USING AN EXPRESSION

mvi = lcomposite.expression(
    '(NIR - GREEN)/(SWIR - GREEN)', {
      'NIR': lcomposite.select('B5_median'),
      'GREEN': lcomposite.select('B3_median'),
      'SWIR': lcomposite.select('B6_median')
}).rename('mvi')

#################################/
###########/   CREATE MVI MASK   ###########/
#################################/

#=================================================================
# STEP 6. THRESHOLDING

# Tweak these values accordingly
# Lower threshold
lower = 4
# Upper threshold
upper = 20

mviina = mvi.lt(upper).add(mvi.gt(lower))
mask = mviina.eq(2)
maskedmvi = mviina.updateMask(mask).rename('mangrove')

mvimaskedcomposite = lcomposite.updateMask(mask)

#################################/
############/   COMPUTE MHI   ############/
#################################/

#=================================================================

#=================================================================
# STEP 8. COMPUTE MHI USING EXPRESSION

#NBR
nbr = mvimaskedcomposite.expression(
    '((NIR - SWIR) / (NIR + SWIR))', {
      'NIR': mvimaskedcomposite.select('B5_median'),
      'SWIR': mvimaskedcomposite.select('B6_median')
}).rename('NBR')
# GCI
gci = mvimaskedcomposite.expression(
    'NIR / GREEN - 1', {
      'NIR': mvimaskedcomposite.select('B5_median'),
      'GREEN': mvimaskedcomposite.select('B3_median')
}).rename('GCI')
# SIPI
sipi = mvimaskedcomposite.expression(
    '((NIR - BLUE) / (NIR - RED))', {
      'NIR': mvimaskedcomposite.select('B5_median'),
      'RED': mvimaskedcomposite.select('B4_median'),
      'BLUE': mvimaskedcomposite.select('B2_median')
}).rename('SIPI')
# ARVI
arvi = mvimaskedcomposite.expression(
    '((NIR - 2*RED + BLUE) / (NIR + 2*RED + BLUE))', {
      'NIR': mvimaskedcomposite.select('B5_median'),
      'RED': mvimaskedcomposite.select('B4_median'),
      'BLUE': mvimaskedcomposite.select('B2_median')
}).rename('ARVI')

# MHI
mhi = mvimaskedcomposite.expression(
    '102.12*NBR - 4.64*GCI + 178.15*SIPI + 159.53*ARVI - 252.39', {
      'NBR': nbr.select('NBR'),
      'GCI': gci.select('GCI'),
      'SIPI': sipi.select('SIPI'),
      'ARVI': arvi.select('ARVI')
}).rename('MHI')

#################################/
#########   PRINTING DATA TO CONSOLE   #########/
#################################/


#################################/
###########   DISPLAY THE DATA   ###########/
#################################/

visParams = {'bands': ['B4_median',  'B3_median',  'B2_median'], 'min': 0, 'max': 0.2}
# Palette for MHI
palette = [
  '000000', 'ff0000', 'ffff00', '00ff00']

# Display L8
Map.addLayer(composite, visParams, 'Landsat 8 Composite')
Map.addLayer(maskedcomposite, visParams, 'composite')

# Display MVI
Map.addLayer(mvi, {}, 'mvi', False)
Map.addLayer(mviina, {}, 'mviina', False)
Map.addLayer(maskedmvi, {}, 'maskedmvi', False)

# Display MHI
Map.addLayer(mhi, {'min': 0, 'max': 100, 'palette': palette}, 'MHI')

#################################/
Map


# In[ ]:




