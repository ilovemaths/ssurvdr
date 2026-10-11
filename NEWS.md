# ssurvdr 0.0.0.9006

* Consolidated the 0.0.0.9004 core and 0.0.0.9005 sparse-path hotfix.
* Declared MASS and survival as required runtime imports; removed unused CVXR and Matrix suggestions.
* Normalised the sole author's name to Kabir Opeyemi Olorede, preserving author and maintainer roles and the existing email.
* Added help pages for the two sparse-path exports and clarified compatibility terminology for the weighting diagnostic.
* Added a short worked vignette and release validation workflow.
* Removed historical console state and research outputs from the source candidate.
* Preserved the numerical estimation and criterion formulas supplied in version 0.0.0.9005.
* Fixed automatic evaluation after screening to map retained predictor indices back to the original matrix. Explicit original-index selection objects used by the simulation driver retain their existing behaviour.
