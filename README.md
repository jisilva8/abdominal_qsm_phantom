# Abdominal QSM Phantom Toolbox

**Abdominal QSM Phantom** is a Julia toolbox for simulating and analyzing abdominal quantitative susceptibility mapping (QSM). The toolbox enables creation of synthetic susceptibility maps, R2* maps, fieldmaps, and simulation of complex water-fat signal acquisitions. It is ideal for validating water-fat separation, background field removal, and QSM methods.

**Main  Features**
   1. Synthetic and editable susceptibility and R2* maps.
   2. Multi-peak fat models.
   3. rBW-dependant chemical shift displacement of fat signal.
   4. rBW-dependant SNR.
   5. Multiecho simulations at user-provided TEs.

---

## Installation and Setup

To prepare the environment and install all required packages, run the provided `setup.jl` script in the project directory. This script activates the project environment and installs all necessary dependencies.

---

## Usage

Example scripts are located in the `examples/` folder. To run an example and generate simulation results:

- **From Julia REPL:**

  1. Start Julia in the project root directory.
  2. Activate the project environment if not already active:
     ```julia
     using Pkg
     Pkg.activate(".")
     ```
  3. Run an example script, for instance:
     ```julia
     include("examples/example04_simulate_acquisition.jl")
     ```

- **From the command line (Bash or terminal):**

  1. Navigate to the project root directory.
  2. Run Julia with the project environment activated and the example included:
     ```
     julia --project=. -e 'include("examples/example04_simulate_acquisition.jl")'
     ```

Running the example scripts will save output files in the corresponding subfolders under the `results/` directory.

---

## Project Structure

The repository contains the following main folders:

```
C:.
├── examples
│ ├── example01_synthetic_chi.jl
│ ├── example02_synthetic_r2star.jl
│ ├── example03_synthetic_fieldmap.jl
│ └── example04_simulate_acquisition.jl
├── results (Created after running examples)
│ ├── ex01_synthetic_chi
│ ├── ex02_synthetic_r2star
│ ├── ex03_synthetic_fieldmap
│ └── ex04_simulate_acquisition
├── src
│ ├── models
│ ├── mri_signal
│ ├── save_data
│ └── tissue_properties
├── Manifest.toml
├── Project.toml
└── setup.jl
```

- `src/`: Core source code including models, MRI signal simulation, data saving, and tissue property modules.
- `examples/`: Example Julia scripts demonstrating various simulations and use cases.
- `results/`: Folder where output files from the example simulations are saved. Running the examples will generate these folders.
- `Project.toml` and `Manifest.toml`: Julia environment files managing dependencies.
- `setup.jl`: Script to initialize the project environment.

---

## Citing

If you have used QSMPhantom toolbox in a scientific publication, we would appreciate citations to the following work:

Silva J, Milovic C, Lambert M, et al. Toward a realistic in silico abdominal phantom for QSM. *Magn Reson Med*. 2023;1-17. doi: 10.1002/mrm.29597

---

## Acknowledgment

This work was funded by Fondecyt 1181057, Fondecyt 1191710, Fondecyt 1231535, Anid/PIA/Anillo ACT192064, and Millennium Institute for Intelligent Healthcare Engineering CN2021_004.  
Dr. Carlos Milovic was supported by Cancer Research UK Multidisciplinary Award C53545/A24348.  
Dr. Cristobal Arrieta was partially funded by ANID Fondecyt Postdoctorado 2019 #3190763.  
Javier Silva was supported by the National Agency for Research and Development ANID/Becas/Doctorado Nacional 21241374.

---

## Bugs and Questions

Please contact Javier Silva Orellana at [jisilva8@uc.cl](mailto:jisilva8@uc.cl)

---

## Disclaimer

THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL JAVIER SILVA OR HIS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, BUSINESS INTERRUPTION; PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; AND LOSS OF USE, DATA OR PROFITS) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
