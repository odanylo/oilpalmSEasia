var step1_s1PreprocessExternalGeometry = function(start, end, roi)
    {
        var collectionDESC =  ee.ImageCollection('COPERNICUS/S1_GRD').filterBounds(roi)
        .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
        .filter(ee.Filter.eq('instrumentMode', 'IW'))
        .select(['VV', 'angle'])
        .filterDate(start, end);
        
        
        var collection2DESC =  ee.ImageCollection('COPERNICUS/S1_GRD').filterBounds(roi)
        .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
        .filter(ee.Filter.eq('instrumentMode', 'IW'))
        .select(['VH', 'angle'])
        .filterDate(start, end);
        
        
        var collectionASC =  ee.ImageCollection('COPERNICUS/S1_GRD').filterBounds(roi)
        .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
        .filter(ee.Filter.eq('instrumentMode', 'IW'))
        .select(['VV', 'angle'])
        .filterDate(start, end);
        
        
        var collection2ASC =  ee.ImageCollection('COPERNICUS/S1_GRD').filterBounds(roi)
        .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
        .filter(ee.Filter.eq('instrumentMode', 'IW'))
        .select(['VH', 'angle'])
        .filterDate(start, end);
        
        // remove edges
        function maskEdge(img) {
          var mask = img.select(0).unitScale(-25, 5).multiply(255).toByte().connectedComponents(ee.Kernel.rectangle(1,1), 250);
          return img.updateMask(mask.select(0));  
        }
        
        // remove stripes
        function stripeMaskVHDESC(im) {
        // VH backscatter of stripes is often less than -30dB
          var mask = ee.Image(0).where(im.select('VH').lte(-16), 1).not();
          return im.updateMask(mask);
        }
        function stripeMaskVVDESC(im) {
        // VH backscatter of stripes is often less than -30dB
          var mask = ee.Image(0).where(im.select('VV').lte(-12), 1).not();
          return im.updateMask(mask);
        }
  
          function stripeMaskVHASC(im) {
        // VH backscatter of stripes is often less than -30dB
          var mask = ee.Image(0).where(im.select('VH').lte(-22), 1).not();
          return im.updateMask(mask);
        }
        function stripeMaskVVASC(im) {
        // VH backscatter of stripes is often less than -30dB
          var mask = ee.Image(0).where(im.select('VV').lte(-15), 1).not();
          return im.updateMask(mask);
        }
        collectionASC = collectionASC.map(stripeMaskVVASC);
        collection2ASC = collection2ASC.map(stripeMaskVHASC);
        
        
        collectionDESC = collectionDESC.map(stripeMaskVVDESC);
        collection2DESC = collection2DESC.map(stripeMaskVHDESC);
        
       var  collection = collectionASC.merge(collectionDESC);
       var  collection2 = collection2ASC.merge(collection2DESC);
        //compute 20th percentile
        var Perc20VV = collection.select('VV').reduce(ee.Reducer.percentile([20]));
        var Perc20VH = collection2.select('VH').reduce(ee.Reducer.percentile([20]));

        //remove images less than perc20
        var MY_filterVV = function(image) {
        var VV_component = image.select('VV');
        var Perc_20_filter = VV_component.gt(Perc20VV);
        return image.updateMask(Perc_20_filter);
        };
        var MY_filterVH = function(image) {
        var VH_component = image.select('VH');
        var Perc_20_filter = VH_component.gt(Perc20VH);
        return image.updateMask(Perc_20_filter);
        };
        
        
        function gamma0(image) {
          var angle = image.select('angle').resample('bicubic');
          return image.select('..').subtract(angle.multiply(Math.PI/180.0).cos().log10().multiply(10.0)).copyProperties(image);
        }
        
        collection = collection.map(MY_filterVV);
        collection2 = collection2.map(MY_filterVH);

        
        var collection_1_1 = collection
        .select("VV");
        var collection_2_1 = collection2
        .select("VH").map(function(im) {return im.addBands(im.int().glcmTexture(3).select(["VH_savg"]))});
       
        var reducer = ee.Reducer.mean();
        var p1VVmean = collection_1_1.reduce(reducer, 4);
        var p1VHmean = collection_2_1.reduce(reducer, 4);
     
        var input = p1VVmean.multiply(10000000).int32()
                            .addBands(p1VHmean.multiply(5000000).int32())
                            .addBands(p1VVsavg.multiply(5000000).int32())
                            .addBands(p1VHsavg.multiply(5000000).int32());
        return input;
 };









var step2_clustering = function()
{
        var savg = ee.ImageCollection("users/OlhaDanylo/global/tropics/s1_savg_2017"),
            s1 = ee.ImageCollection("users/OlhaDanylo/global/tropics/s1_2017"),
            grid = ee.FeatureCollection("users/OlhaDanylo/global/tropics/tropics_grid"),
            sumatra = /* color: #d63000 */ee.Geometry.Point([104.57957198411214, -3.5075359534240467]),
            dem = ee.Image("USGS/SRTMGL1_003"),
            kalimantan = /* color: #d63000 */ee.Geometry.Point([110.80290172787522, -1.0002196683841287]),
            sumatra_north = /* color: #98ff00 */ee.Geometry.Point([101.10662231547303, 0.7707982892375213]),
            sumatra_far_north = /* color: #0b4a8b */ee.Geometry.Point([98.36004028422303, 2.857039146683248]),
            kalim_north = /* color: #ffc82d */ee.Geometry.Point([112.79607544047303, 1.8252002632183957]),
            kalim_north_west = /* color: #00ffff */ee.Geometry.Point([116.22380981547303, 2.39609937533162]),
            kalim_south_west = /* color: #bf04c2 */ee.Geometry.Point([115.67449340922303, -1.7557590754775945]),
            malaysia = /* color: #ff0000 */ee.Geometry.Point([102.51287231547303, 3.5809878866227494]),
            north_sulawesi = /* color: #98ff00 */ee.Geometry.Point([121.7880285843579, 0.7982014067499168]),
            kimbe = /* color: #ff0000 */ee.Geometry.Point([150.7008338585705, -5.411313694948916]),
            thailand = /* color: #0b4a8b */ee.Geometry.Point([98.02724016577213, 7.714441828067812]),
            north_malaysia = /* color: #ffc82d */ee.Geometry.Point([101.51723461406596, 5.82883444168634]),
            far_north_kalimantan = /* color: #00ffff */ee.Geometry.Point([116.90294520580426, 5.695498763356799]),
            papua = /* color: #bf04c2 */ee.Geometry.Point([142.46581877974631, -3.889439202525224]),
            danylo_v2 = ee.ImageCollection("users/OlhaDanylo/palmoil/danylo_v2");
            
    
        var ress = ee.Image(0);
       
        s1 = s1.mosaic();
        var input = s1.select("VH").addBands(s1.select("VV").multiply(1.25))
                          .addBands(s1.select("VH").subtract(s1.select("VV")).multiply(10))
                          .addBands(s1savg.mosaic());
        
        
        ////// ###################
        ////// ################### south sumatra
        ////// ###################
        
        var region = grid.filterBounds(sumatra);
        
        var bands = input.bandNames();
        
        var training = input.sample({
            region: region,
            scale: 20,
            numPixels: 50000
        });
        
        var clusterers = [
          ee.Clusterer.wekaXMeans(10, 15)
        ];
        
        var clus2 = clusterers[0].train(training);
        var res = input.cluster(clus2);
        
        ress = ress.add(res.eq(6).or(res.eq(5)));
        
        
        ////// ###################
        ////// ################### kalimantan
        ////// ###################
        
        region = grid.filterBounds(kalimantan);
        bands = input.bandNames();
        
        training = input.sample({
            region: region,
            scale: 20,
            numPixels: 50000
        });
        
        clusterers = [
          ee.Clusterer.wekaXMeans(10, 15)
        ];
        
        clus = clusterers[0].train(training);
        res = input.cluster(clus);
        
        ress = ress.add(res.eq(3).or(res.eq(5)).or(res.eq(4)).or(res.eq(4)));
        
        
        
        ////// ###################
        ////// ################### malaysia
        ////// ###################
        
        region = grid.filterBounds(malaysia);
        bands = input.bandNames();
        
        
        training = input.sample({
            region: region,
            scale: 20,
            numPixels: 50000,
            seed: 11
        });
        
        clus = clusterers[0].train(training);
        res = input.cluster(clus);
        
        ress = ress.add(res.eq(9)
        .or(res.eq(10))
        .or(res.eq(3))
        .or(res.eq(4)));
        
        
        
        ////// ###################
        ////// ################### north west kalimantan
        ////// ###################
        
        region = grid.filterBounds(kalim_north_west);
        
        bands = input.bandNames();
        
        training = input.sample({
            region: region,
            scale: 20,
            numPixels: 50000,
            seed: 11
        });
        
        
        clus = clusterers[0].train(training);
        res = input.cluster(clus);
        
        ress = ress.add(res.eq(13));
        
        
        
        ////// ###################
        ////// ################### kalimantan north
        ////// ###################
        
        region = grid.filterBounds(kalim_north);
        
        bands = input.bandNames();
        
        
        training = input.sample({
            region: region,
            scale: 20,
            numPixels: 50000,
            seed: 11
        });
        
        
        clus = clusterers[0].train(training);
        res = input.cluster(clus);
        
        ress = ress.add(res.eq(10).or(res.eq(9)));
        
        
        
        ///// ###################
        ////// ################### sumatra far north
        ////// ###################
        
        region = grid.filterBounds(sumatra_far_north);
        bands = input.bandNames();
        
        
        training = input.sample({
            region: region,
            scale: 20,
            numPixels: 50000,
            seed: 11
        });
        
        clus = clusterers[0].train(training);
        res = input.cluster(clus);
        
        ress = ress.add(res.eq(11).or(res.eq(0)));
        
        
        
        ////// ###################
        ////// ################### kalimantan south west
        ////// ###################
        
        region = grid.filterBounds(kalim_south_west);
        bands = input.bandNames();
        
        training = input.sample({
            region: region,
            scale: 20,
            numPixels: 50000
        });
        
        
        clus = clusterers[0].train(training);
        res = input.cluster(clus);
        
        ress = ress.add(res.eq(3));
      
        
        ////// ###################
        ////// ################### thailand
        ////// ###################
        
        region = grid.filterBounds(thailand);
        
        bands = input.bandNames();
        
        training = input.sample({
            region: region,
            scale: 20,
            numPixels: 50000
        });
        
        clus = clusterers[0].train(training);
        res = input.cluster(clus);
        
        ress = ress.add(res.eq(8));
       
        ////// ###################
        ////// ################### north malaysia
        ////// ###################
        
        region = grid.filterBounds(north_malaysia);

        // Use these bands in the prediction.
        bands = input.bandNames();
        
        
        training = input.sample({
            region: region,
            scale: 20,
            numPixels: 50000
        })
        
        var clus = clusterers[0].train(training);
        res = input.cluster(clus);
        
        ress = ress.add(res.eq(3));
       
        ////// ###################
        ////// ################### far north kalimantan
        ////// ###################
        
        region = grid.filterBounds(far_north_kalimantan);
        bands = input.bandNames();
        
        
        training = input.sample({
            region: region,
            scale: 20,
            numPixels: 50000
        });
        
        clus = clusterers[0].train(training);
        res = input.cluster(clus);
        
        ress = ress.add(res.eq(3).or(res.eq(11)));
        
        
        ////// ###################
        ////// ################### north sulawesi
        ////// ###################
        
        region = grid.filterBounds(north_sulawesi);
        bands = input.bandNames();
        
        training = input.sample({
            region: region,
            scale: 20,
            numPixels: 50000
        });
        
        clus = clusterers[0].train(training);
        res = input.cluster(clus);
        
        
        
        ////// ###################
        ////// ################### papua
        ////// ###################
        
        region = grid.filterBounds(papua);

        bands = input.bandNames();
        
        
        training = input.sample({
            region: region,
            scale: 20,
            numPixels: 50000
        });
        
        clus = clusterers[0].train(training);
        res = input.cluster(clus);
        
        ress = ress.add(res.eq(8));
        
        
        ////// ###################
        ////// ################### papua
        ////// ###################
        
        region = grid.filterBounds(kimbe);
        bands = input.bandNames();
        
        training = input.sample({
            region: region,
            scale: 20,
            numPixels: 50000
        });
        
        
        clus = clusterers[0].train(training);
        res = input.cluster(clus);
        
        ress = ress.add(res.eq(6));
        
        // Greenest Pixel Composite
        function addNDVI(input) {
          var ndvi = input.normalizedDifference(['B5', 'B4']).rename("ndvi");
          return input.addBands(ndvi);
        }
        
        var collection = ee.ImageCollection("LANDSAT/LC08/C01/T1_TOA")
            .filterDate('2017-01-01', '2018-12-01');
        
        var ndviCollection = collection.map(addNDVI);
        var composite = ndviCollection.qualityMosaic('ndvi');
        
        var palmoilClassification = ress.mask(ress.gte(7).and(composite.select("ndvi").gte(0.5)).and(ee.Terrain.slope(dem).lte(10)))
                                .gte(7);
        
        Export.image.toAsset(palmoilClassification, "danylo");
};






