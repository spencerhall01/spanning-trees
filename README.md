# Spanning Trees on Quantum Annealers
   
   This repository contains MATLAB scripts used in studies of random spanning trees on quantum annealers with a Chimera topology (*e.g.*, old D-Wave systems). These scripts were developed as part of collaborative research efforts focused on performance evaluation of quantum hardware using planted solutions and mirror symmetry.
   
   Core scripts include random spanning tree generation, Chimera graph modeling, Ising problem formulation, and interaction/field assignment. These scripts are old and are probably no longer fully functional on current machines; they are presented here for archival and educational purposes. The core ideas and algorithms may still be useful.
   
## Repository Structure
   
    .
    ├── archive/     # Scripts as they were at the end of the research projects
    ├── core/        # Core scripts used during the studies - may be updated later
    ├── data/        # I/O files for simple examples or testing - empty for now
    ├── docs/        # Notes, images, TODOs, links, etc. - empty for now
    ├── LICENSE      # License file (MIT)
    └── README.md    # Project overview and documentation
   
## Contents
   
   ### `core/` – Core research scripts
   
   These are the most complete and reusable scripts, designed for Ising model studies on early D-Wave machines. Everything here has been tested and was in a working state at one point. If you would like to adapt the code, this would be the place to start. Copies of these are in the archive; the copies here may be updated in the future.
   
   - `randspan.m` – Constructs a spanning tree via random walk on a subgraph of the D-Wave Chimera architecture and plants the Ising problem. The spanning tree algorithm is based on the method introduced by Broder (1989); see [related publications](#related-research-publications) below.
      
   - `randspanmir.m` – Same as above, but enforces vertical mirror symmetry as an additional constraint on the Ising problem. This technique was adapted and extended in follow-up work by other collaborators; see [related publications](#related-research-publications) below.
   
   **NOTE:** These scripts are monolithic; this was before I really learned about good software practices like modularization. However, I think the code is documented well enough that one should be able to follow the process. Some parts are commented out because I wasn't using them, not because they weren't functional.
   
   ### `archive/` – Scripts as they were when the research stopped
		
   Includes eleven MATLAB scripts reflecting different stages of the development and research process. The core scripts incorporate most of the functionality of these scripts. Nothing in these files should be modified.
   
   **WARNING:** These scripts are unfactored and contain some incomplete or incorrect implementations. They are preserved here as a historical snapshot and for later reference.
   
   ### `data/` - Small data files
   
   This folder is for storing *small* data files used in simple examples or during testing. It's also a good place to keep input parameters, particular planted solutions, *etc.*, if you want to load them from a file. You may also use this folder to store your own output locally, but **large data files should NOT be uploaded to the repository**.
   
   ### `docs/` - Related documents

   This folder is for storing simple documents like notes, TODOs, small images, links to resources, *etc.*, related to the project. Again, **massive files should NOT be uploaded to the repository**.
   
## Related research publications
   
   These provide background and context for the code, including theoretical motivations and performance analyses.
   
   - Broder. "Generating random spanning trees." [doi:10.1109/SFCS.1989.63516](https://doi.org/10.1109/SFCS.1989.63516)
   - Hall *et al.* "A Study of Spanning Trees on a D-Wave Quantum Computer." [doi:10.1016/j.phpro.2015.07.109](https://doi.org/10.1016/j.phpro.2015.07.109)
   - Hobl. "Simulating on the D-Wave Two and emulating its behavior on an ordinary computer." [https://juser.fz-juelich.de/record/276310](https://juser.fz-juelich.de/record/276310)
   - Novotny, Hobl, Hall, and Michielsen. "Spanning Tree Calculations on D-Wave 2 Machines." [doi:10.1088/1742-6596/681/1/012005](https://doi.org/10.1088/1742-6596/681/1/012005)
   - Hall. "A Study Of The Performance Of D-Wave Quantum Computers Using Spanning Trees." [https://scholarsjunction.msstate.edu/td/302/](https://scholarsjunction.msstate.edu/td/302/)
   - Perera, Hall, and Novotny. “Validating the solutions of the D-Wave quantum annealers through graph mirroring.” [http://meetings.aps.org/link/BAPS.2016.MAR.F45.4](http://meetings.aps.org/link/BAPS.2016.MAR.F45.4)
   - Perera and Novotny. "An answer checking method for quantum annealers." [doi:doi.org/10.1088/1742-6596/750/1/012005](https://doi.org/10.1088/1742-6596/750/1/012005)
   - Perera *et al.* "Benchmarking quantum annealers using symmetries in embedded subgraphs." [doi:10.48550/arXiv.1708.01026](https://doi.org/10.48550/arXiv.1708.01026)
   - Bhalgamiya *et al*. "Quantum Annealing Error Mitigation Using Mirror Symmetries on Different Generations of Quantum Annealers." [doi:10.1109/QCE53715.2022.00140](https://doi.org/10.1109/QCE53715.2022.00140)
   
## Requirements
   
   - MATLAB (tested on older versions)
   - Access to D-Wave hardware (or simulator)
   - Valid API token for remote execution
   - Basic familiarity with Ising models and quantum annealers
   
## Usage
   
   1. Edit the scripts to set:
      - Solver name (e.g., `C12`)
      - API credentials
      - Chip layout (*x* × *y* unit cells)
      - Rectangular subgraph size and position
      
   2. Run `randspan.m` or `randspanmir.m` in MATLAB.
   
   3. Review output files like `J.dat`, `average_magnetization.dat`, and any exported energy/magnetization statistics.
   
   See comments in each file for details. Many hardcoded parameters can be made interactive or *vice versa*. Some parts are commented out because I wasn't using them at the time the research ended.
   
## License
   
   This repository is licensed under the MIT License. See [`LICENSE`](./LICENSE) for details.
   
## Citation
   
   If you use this code or portions of it, please cite or acknowledge this repository :)
   