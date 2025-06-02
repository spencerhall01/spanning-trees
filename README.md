# Spanning Trees on Quantum Annealers
   
   This repository contains MATLAB scripts used in studies of random spanning trees on quantum annealers with a Chimera topology (*e.g.*, old D-Wave systems). These scripts were developed as part of collaborative research efforts focused on performance evaluation of quantum hardware using planted solutions and mirror symmetry.
   
   Core scripts include random spanning tree generation, Chimera graph modeling, Ising problem formulation, and interaction/field assignment. These scripts are old and probably no longer fully functional on current machines; they are presented here for archival and educational purposes. The core ideas and algorithms may still be useful.
   
## Repository Structure
   
    .
    ├── archive/     # Scripts as they were at the end of the research projects
    ├── core/        # Core scripts used during the studies - may be updated later
    ├── data/        # I/O files for simple examples or testing - empty for now
    ├── docs/        # Notes, images, TODOs, links, etc. - empty for now
    ├── LICENSE      # License file (MIT)
    └── README.md    # Project overview and documentation
   
## Contents
   
   ### core/ – Core research scripts
   
   These are the most complete and reusable scripts, designed for Ising model studies on early D-Wave machines. Everything here has been tested and was in a working state at one point. If you would like to adapt the code, this would be the place to start. Copies of these are in the archive; the copies here may be updated in the future.
   
   * `randspan.m` – Constructs a spanning tree via random walk on a subgraph of the D-Wave Chimera architecture and plants the Ising problem. Formats and saves the output.
      
   * `randspanmir.m` – Same as above, but with vertical mirror symmetry added as an additional constraint. Horizontal symmetry was implemented later in a separate script by another student researcher (see Hobl in `./docs`).
   
   **NOTE:** These scripts are monolithic; this was before I really learned about good software practices like modularization. However, I think the code is documented well enough that one should be able to follow the process. Some parts are commented out because I wasn't using them, not because they weren't functional.
   
   ### archive/ – Scripts as they were when the research stopped
		
   Includes ten MATLAB scripts reflecting different stages of the development and research process. The core scripts incorporate most of the functionality of these scripts. Nothing in these files should be modified.
   
   **WARNING:** These scripts are unfactored and contain some incomplete or incorrect implementations. They are preserved here as a historical snapshot and for later reference.
   
   ### data/ - Small data files
   
   This folder is for storing *small* data files used in simple examples or during testing. It's also a good place to keep input parameters, particular planted solutions, *etc.* if you want to load them from a file. You may also use this folder to store your own output locally, but **large data files should NOT be uploaded to the repository**.
   
   ### docs/ - Related documents

   This folder is for storing simple documents like notes, TODOs, small images, links to resources, *etc.* related to the project. Again, **massive files should NOT be uploaded to the repository**.
   
## Related research publications
   
   These provide background and context for the code, including theoretical motivations and performance analyses.
   
   - [Hall *et al.* "A Study of Spanning Trees on a D-Wave Quantum Computer"](https://doi.org/10.1016/j.phpro.2015.07.109)  
   - [Novotny *et al.* "Spanning Tree Calculations on D-Wave 2 Machines"](https://doi.org/10.1088/1742-6596/681/1/012005)  
   - [Hobl. "Evaluation of Spanning Tree Sampling on a Quantum Annealer" (Master’s Thesis)](https://juser.fz-juelich.de/record/276310)
   
## Requirements
   
   - MATLAB (tested on R2008b)
   - Access to D-Wave hardware (or simulator)
   - Valid API token for remote execution
   - Basic familiarity with Ising models and quantum annealers
   
## Usage
   
   ### How to Run the Scripts
   
   1. Edit the scripts to set:
      * Solver name (e.g., `C12`)
      * API credentials
      * Chip layout (*x* × *y* unit cells)
      * Rectangular subgraph size and position
      
   2. Run `randspan.m` or `randspanmir.m` in MATLAB.
   
   3. Review output files like `J.dat`, `average_magnetization.dat`, and any exported energy/magnetization statistics.
   
   See comments in each file for details. Many hardcoded parameters can be made interactive or *vice versa*. Some parts are commented out because I wasn't using them at the time the research ended.
   
## License
   
   This repository is licensed under the MIT License. See [`LICENSE`](./LICENSE) for details.
   
## Citation
   
   If you use this code or portions of it, please cite or acknowledge :)
	