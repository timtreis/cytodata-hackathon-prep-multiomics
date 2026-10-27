import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import spatialdata as sd
from spatialdata.models import TableModel

import prep


def _sdata() -> sd.SpatialData:
    obs = pd.DataFrame(
        {
            "Well": ["A1", "A1", "A2", "A2"],
            "WellLabel": ["Cellecta", "Cellecta", "CellectaTNFa", "CellectaTNFa"],
            "cell_id_global": [1, 2, 3, 4],
            "region_key": pd.Categorical(["A1-labels-Cell", "A1-labels-Cell", "A2-labels-Cell", "A2-labels-Cell"]),
            "Cell": [1, 2, 1, 2],
        },
        index=list("abcd"),
    )
    rna = sp.csr_matrix(np.array([[300, 0, 10], [50, 0, 0], [400, 5, 1], [500, 0, 2]]))
    img = np.ones((4, 3), dtype="float32")
    img[1, 2] = np.nan  # only a _NuclearRatio NaN: must not flag
    img[3, 0] = np.nan  # real NaN: must flag
    geom = pd.DataFrame({"NuclearArea": [10, 10, 0, 10]}, index=obs.index)
    art = pd.DataFrame({"any": [0, 0, 0, 1]}, index=obs.index)
    ta = pd.DataFrame(
        {"MostAbundantTarget": ["TGFBR2_1", "NT_Human_CRko_BB_2", "OR1I1_3", "Undetermined"], "Confidence": [0.9, 0.6, 0.4, 0.0]},
        index=obs.index,
    )

    def table(X, var, **obsm):
        a = ad.AnnData(X, obs=obs.copy(), var=pd.DataFrame(index=var), obsm={"cell_geometry": geom, "spatial_artifacts": art, **obsm})
        return TableModel.parse(a, region=["A1-labels-Cell", "A2-labels-Cell"], region_key="region_key", instance_key="Cell")

    meta = ad.AnnData(uns={"wells": pd.DataFrame({"cell_type": ["HELA", "HELA"]}, index=["A1", "A2"])})
    return sd.SpatialData(
        tables={
            "global-table-ThreePrimeUntargeted": table(rna, ["GAPDH", "MT-CO1", "ACTB"]),
            "global-table-imaging": table(img, ["Actin", "Nucleus", "Actin_NuclearRatio"]),
            "global-table-B01-Cellecta": table(sp.csr_matrix((4, 1)), ["g"], target_assignment=ta),
            "global-table-B01-QC_Panel": table(sp.csr_matrix((4, 1)), ["q"]),
            "global-table-metadata": meta,
        }
    )


def test_prep():
    s = _sdata()
    ann = prep.annotations(s)
    assert ann["condition"].tolist() == ["baseline", "baseline", "TNFa", "TNFa"]
    assert ann["target_gene"].tolist() == ["TGFBR2", "NTC", "OR1I1", "Undetermined"]
    assert ann["guide_assigned"].tolist() == [True, True, False, False]
    assert ann["qc_low_rna"].tolist() == [False, True, False, False]
    assert ann["qc_imaging_nan"].tolist() == [False, False, False, True]
    assert ann["qc_no_nucleus"].tolist() == [False, False, True, False]
    assert ann["pass_qc"].tolist() == [True, False, False, False]

    prep.annotate(s, ann)
    wells = s.tables["global-table-metadata"].uns["wells"]
    assert (wells["cell_type"] == "A549").all() and (wells["cell_type_original"] == "HELA").all()

    sub = prep.subset(s, ["A2"])
    t = sub.tables["global-table-imaging"]
    assert t.obs["Well"].tolist() == ["A2", "A2"] and t.uns["spatialdata_attrs"]["region"] == ["A2-labels-Cell"]
