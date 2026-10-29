import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from model_he_holdout import evaluate

def make(has_history,seed):
    rng=np.random.default_rng(seed);groups={}
    for i in range(90):
        current=rng.normal(size=(8,4));past=current+rng.normal(size=(8,4))
        loc=rng.normal(0,10,size=(8,2));n=rng.integers(15,45,size=(8,1))
        x=np.column_stack([past,current,loc,n])
        hist=np.linalg.norm(past-current,axis=1)
        y=.25+.07*current[:,0]+rng.normal(0,.01,8)
        if has_history:y=y+.07*hist
        groups[str(i)]={"species":str(i),"role":i%2,"x":x,"y":y}
    return groups
class TestSyntheticKnownTruth(unittest.TestCase):
    def test_added_history_improves_when_true(self):
        r=evaluate(make(True,1234),draws=99)
        self.assertGreater(r["relative_heldout_MSE_improvement"],.5)
        self.assertTrue(r["positive_prediction_gain_95lower_gt_zero"])
    def test_null_history_does_not_appear_positive(self):
        r=evaluate(make(False,5421),draws=99)
        self.assertFalse(r["positive_prediction_gain_95lower_gt_zero"])
if __name__=="__main__":unittest.main()
