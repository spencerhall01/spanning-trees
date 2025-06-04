# Spanning Trees on Quantum Annealers
   
   This repository contains MATLAB scripts used in studies of random spanning trees on quantum annealers with a Chimera topology (_e.g._, old D-Wave systems). These scripts were developed as part of collaborative research efforts focused on performance evaluation of quantum hardware using planted solutions and mirror symmetry.
   
   Core scripts include random spanning tree generation, Chimera graph modeling, Ising problem formulation, and interaction/field assignment. These scripts are old and are probably no longer fully functional on current machines; they are presented here for archival and educational purposes. The core ideas and algorithms may still be useful.
   
## Repository Structure
   
    .
    ├── archive/     # Scripts as they were at the end of the research projects
    ├── data/        # I/O files for simple examples or testing - empty for now
    ├── docs/        # Notes, images, TODOs, links, etc. - empty for now
    ├── scripts/     # Core scripts used during the studies - may be updated later
    ├── LICENSE      # License file (MIT)
    └── README.md    # Project overview and documentation
   
## Contents
   
   ### `scripts/` – Core research scripts
   
   These are the most complete and reusable scripts, designed for Ising model studies on early D-Wave machines. Everything here has been tested and was in a working state at one point. If you would like to adapt the code, this would be the place to start. Copies of these are in the archive; the copies here may be updated in the future.
   
   - `randspan.m` – Constructs a spanning tree via random walk on a subgraph of the D-Wave Chimera architecture and plants the Ising problem. The spanning tree algorithm is based on the method introduced by Broder (1989); see [Related Publications](#related-publications) below.
      
   - `randspanmir.m` – Same as above, but enforces vertical mirror symmetry as an additional constraint on the Ising problem. This technique was adapted and extended in follow-up work by other collaborators; see [Related Publications](#related-publications) below.
   
   **NOTE:** These scripts are monolithic; this was before I really learned about good software practices like modularization. However, I think the code is documented well enough that one should be able to follow the process. Some parts are commented out because I wasn't using them, not because they weren't functional.
   
   ### `archive/` – Scripts as they were when the research stopped
		
   Includes eleven MATLAB scripts reflecting different stages of the development and research process. The core scripts incorporate most of the functionality of these scripts. Nothing in these files should be modified.
   
   **WARNING:** These scripts are unfactored and contain some incomplete or incorrect implementations. They are preserved here as a historical snapshot and for later reference.
   
   ### `data/` - Small data files
   
   This folder is for storing *small* data files used in simple examples or during testing. It's also a good place to keep input parameters, particular planted solutions, _etc._, if you want to load them from a file. You may also use this folder to store your own output locally, but **large data files should NOT be uploaded to the repository**.
   
   ### `docs/` - Small documents and images

   This folder is for storing simple documents like notes, TODOs, small images, resources, _etc._, related to the project. Again, **massive files should NOT be uploaded to the repository**.
   
## Related Publications
   
   These provide background and context for the code, including theoretical motivations and performance analyses.
   
   ### Publications related to spanning trees
   - Broder. "Generating random spanning trees." _30th Annual Symposium on Foundations of Computer Science_, 1989. [doi:10.1109/SFCS.1989.63516](https://doi.org/10.1109/SFCS.1989.63516)
   - **Hall** _et al._ "A Study of Spanning Trees on a D-Wave Quantum Computer." _Physics Procedia_, 2015. [doi:10.1016/j.phpro.2015.07.109](https://doi.org/10.1016/j.phpro.2015.07.109)
   - Hobl. "Simulating on the D-Wave Two and emulating its behavior on an ordinary computer." Master's thesis, 2015. [https://juser.fz-juelich.de/record/276310](https://juser.fz-juelich.de/record/276310)
   - Novotny, Hobl, **Hall**, and Michielsen. "Spanning Tree Calculations on D-Wave 2 Machines." _Journal of Physics: Conference Series_, 2016. [doi:10.1088/1742-6596/681/1/012005](https://doi.org/10.1088/1742-6596/681/1/012005)
   - **Hall**. "A Study Of The Performance Of D-Wave Quantum Computers Using Spanning Trees." Master's thesis, 2018. [https://scholarsjunction.msstate.edu/td/302/](https://scholarsjunction.msstate.edu/td/302/)

   ### Publications related to mirror symmetry
   - Perera, **Hall**, and Novotny. “Validating the solutions of the D-Wave quantum annealers through graph mirroring.” APS March Meeting, 2016. [http://meetings.aps.org/link/BAPS.2016.MAR.F45.4](http://meetings.aps.org/link/BAPS.2016.MAR.F45.4)
   - Perera and Novotny. "An answer checking method for quantum annealers." _Journal of Physics: Conference Series_, 2016. [doi:doi.org/10.1088/1742-6596/750/1/012005](https://doi.org/10.1088/1742-6596/750/1/012005)
   - Perera _et al._ "Benchmarking quantum annealers using symmetries in embedded subgraphs." arXiv, 2022. [doi:10.48550/arXiv.1708.01026](https://doi.org/10.48550/arXiv.1708.01026)
   - Bhalgamiya _et al._ "Quantum Annealing Error Mitigation Using Mirror Symmetries on Different Generations of Quantum Annealers." _IEEE International Conference on Quantum Computing and Engineering (QCE)_, 2022. [doi:10.1109/QCE53715.2022.00140](https://doi.org/10.1109/QCE53715.2022.00140)
   
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
   