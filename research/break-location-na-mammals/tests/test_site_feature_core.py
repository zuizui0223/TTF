import importlib,importlib.util,unittest,tempfile,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import LineString,box

class Tests(unittest.TestCase):
    def module(self):
        spec=importlib.util.find_spec('site_feature_core')
        self.assertIsNotNone(spec,'site_feature_core extraction functions are not implemented')
        return importlib.import_module('site_feature_core')
    def test_antimeridian_is_short_not_zero_meridian(self):
        m=self.module();v=m.dense_vertices([LineString([(179.9,30),(-179.9,30)])])
        self.assertTrue(np.all(np.abs(v[:,0])>179.8))
        self.assertGreater(len(v),10)
    def test_preserve_internal_ecotone(self):
        m=self.module();v=m.dense_vertices([box(0,0,1,1),box(1,0,2,1)])
        d=m.vertex_distances([(1,.5),(.5,.5)],v)
        self.assertLess(d[0],.8);self.assertGreater(d[1],50.)
    def test_sphere_great_circle_distance(self):
        m=self.module();d=m.vertex_distances([(0,0)],np.array([[1.,0.]]))
        self.assertAlmostEqual(d[0],111.19508,places=3)
    def test_masked_raster_and_oob_are_nan(self):
        m=self.module()
        with tempfile.TemporaryDirectory() as p:
            path=Path(p)/'x.tif'
            with rasterio.open(path,'w',driver='GTiff',height=2,width=2,count=1,dtype='float32',crs='EPSG:4326',transform=from_origin(0,2,1,1),nodata=-9999) as ds:
                ds.write(np.array([[0.,-9999],[-3.,7.]],'f4'),1)
            v=m.raster_values(path,[(.5,1.5),(1.5,1.5),(.5,.5),(3.,3.)])
        np.testing.assert_equal(v,np.array([0.,np.nan,-3.,np.nan]))
    def test_feature_missingness_is_not_imputed(self):
        m=self.module();x=np.ones((12,13));x[1,0]=np.nan
        self.assertEqual(m.complete_rows(x).tolist(),[True,False]+[True]*10)
    def test_bad_coords_raise(self):
        m=self.module()
        with self.assertRaises(ValueError):m.sphere_xyz([[181,0]])
if __name__=='__main__':unittest.main()
