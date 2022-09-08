# --- VERSIONS OF CSDL AND CSDL_OM FOR lsdo_mesh & motor_project: ---
- use git checkout <commit-id> to change to these versions in DETACHED_HEAD state (recommended)
    - git checkout <commit-id> . will leave you with all future commits (not recommended)
- use git switch - to revert back to main master branch
- csdl: *https://github.com/LSDOlab/csdl/tree/98d7f6d43a74042a5f8ded2a2b20d431f0c1c342*
    - commit-id: 98d7f6d
- csdl_om: *https://github.com/LSDOlab/csdl_om/tree/80e67db290493054d172a4813d3f16c28bbdca50*
    - commit-id: 80e67db

# --- TO-DO ---
- fix parametrization of points for rotated instances (appearance of NaNs)
    - the issue currently in the FFD faces appears for a delta applied to theta, where the FFD face will cross over \theta=\pi, and this messes up the FFD face parametrization (run the test_magnet_FFD.py case and zoom into the magnet just above the "-" x-axis, and the theta delta is improperly computed)
    - **UPDATE**: FFD face parametrization across "-" x-axis is fixed; currently working on edges
    - **UPDATE (09/07)**: edge parametrization works now across "-" x-axis for magnet test
- code clean-up
- generalization for cartesian (x-y) coordinates as well
    - implementation is only in polar right now

- interface current magnet FFD test with Ru's code
    - just to check the FFD-related stuff

- expand the FFD faces to consider other parts of the motor mesh and put into a shape parameter model


- COME UP WITH A CONSISTENT NAMING SCHEME FOR EDGES AND VERTICES TO AVOID CONFUSION
    - mesh vertices: user-defined points (making up major subdomains like the 4 points defining a magnet)
    - edge nodes: nodes along edges between vertices
    - mesh nodes: internal mesh nodes (deformation done by Ru's FEniCS mesh manipulation code)

# --- NOTES ---
