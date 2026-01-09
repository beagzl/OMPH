# Optimal Monetary Policy with Heterogeneity

A MATLAB/Dynare toolkit to solve continuous-time heterogeneous-agent models in the sequence space and to compute optimal policies within those models, following the methodology of the paper "Firm Heterogeneity, Capital Misallocation and Optimal Monetary Policy”, by González, Nuño and Thaler (2026)

Copyright (c) 2026 Beatriz González, Galo Nuño, Dominik Thaler

If you use this code in your research, please cite:

```bibtex
@article{gonzalez2026firm,
  title={Firm Heterogeneity, Capital Misallocation and Optimal Monetary Policy},
  author={Gonz{\'a}lez, Beatriz and Nu{\~n}o, Galo and Thaler, Dominik},
  year={2026}
}
```


## 1- Overview

This repository provides automated code generation for solving HANK and RANK models using Dynare. The toolkit:
- **Automatically derives** first-order conditions for Ramsey optimal policy
- **Generates discretized** Dynare model files from continuous-time equations
- **Handles wealth distribution** dynamics (Kolmogorov Forward Equation for HANK)
- **Supports multiple policy** regimes: Ramsey optimal, Taylor rule, fixed inflation, exogenous paths
- **Analyzes transitional dynamics** with various shock types (monetary, preference, cost-push, TFP, credit)
 
## 2- Requirements

- MATLAB R2019b or later
- Symbolic Math Toolbox
- Dynare 4.6 or later
  - Ensure Dynare is installed and on your MATLAB path

## 3- Usage

### 3.1- Features

**Model Types**
- **HANK**: Full heterogeneous agent model with wealth distribution dynamics
- **RANK**: Representative agent benchmark for comparison

**Policy Regimes** (controlled via `OMP` parameter)
- `OMP=0`: Taylor rule with interest rate smoothing
- `OMP=1`: Ramsey optimal monetary policy
- `OMP=2`: Fixed inflation target (π = π̄)

**Shock Types**
- Monetary policy shocks
- Time preference (discount rate) shocks
- Cost-push (markup) shocks
- TFP (productivity) shocks

**Solution Methods**
- Linear approximation
- Nonlinear perfect foresight

### 3.2- Quick Start

 **Create the folder "/_simulations"** inside the main folder
 
 **Run the Master file** `MASTER.m`

1. **Generate both models** (run just once):
   ```matlab
   % Run Section 1 of MASTER.m
   % This creates RANK and HANK model files
   ```

2. **Configure simulation** in `model_parameters.m`:
   ```matlab
   OMP = 1;              % Ramsey optimal policy
   shocktype = 2;        % Time preference shock
   shocksize = -100;     % -100 basis points
   ```

3. **Run simulations**:
   ```matlab
   % Run Section 2 of MASTER.m
   % Results automatically saved to .mat files in "/_simulations" folder
   ```

### 3.3- Detailed Workflow

#### Step 1: Model Generation

The `buildfile.m` script generates Dynare code from symbolic equations:

```matlab
modelfilename = 'HANK';  % Choose model type
buildfile                % Generate Dynare files
```

**Generated files:**
- `HANK_DynamicsEquCond.mod` - Private equilibrium dynamics
- `HANK_DynamicsRamsey.mod` - Ramsey FOCs
- `HANK_SSmultfile.m` - Steady-state multiplier calculations
- `HANK_DeclareVars.mod` - Variable declarations
- `HANK_DeclareMult.mod` - Multiplier declarations
- `HANK_SSfilecall.mod` - Steady-state calculations

#### Step 2: Configure Simulation

Edit `model_parameters.m` to set:

**Model selection:**
```matlab
modelfilename = 'HANK';  % or 'RANK'
N = 100;                 % Wealth grid size (HANK only)
```

**Policy regime (`OMP`):**
- `0` - Taylor rule
- `1` - Ramsey optimal policy
- `2` - Fixed inflation: `π = π̄`


**Shock configuration:**
```matlab
shocktype = 2;      % 0=none, 1=monetary, 2=preference, 3=cost-push, 4=TFP
shocksize = -100;    % Size in basis points (100 bp = 1%)
```

**Solution method:**
```matlab
linearize = 0;      % 0=nonlinear, 1=linear, 2=linear then nonlinear
```

**Simulation horizon:**
```matlab
N_period = 600;     % Number of periods
Deltat = 1/12;      % Time step (months)
```

#### Step 3: Run Simulation

```matlab
run_mod
```

This script:
1. Loads parameters from `model_parameters.m`
2. Calls Dynare with appropriate macros
3. Solves for steady state
4. Computes transitional dynamics
5. Generates impulse response figures



### Example: Time Preference Shock under Optimal Policy

```matlab
% Configure simulation
modelfilename = 'HANK';
N = 100;
OMP = 1;              % Ramsey optimal policy
shocktype = 2;        % Time preference shock
shocksize = -100;      % -100 basis points


% Generate model (only needed once)
buildfile

% Run Master file
run_mod
```

This simulates a -100bp shock to the time preference rate under optimal monetary policy.


## 4- Customization

### Modifying the Model

To modify the model:

1. Edit equations in `buildfile.m` in the appropriate `case` block and adjust number of forward/backlooking dynamic equations. The user needs to provide the dynamic private equilibrium conditions and the planner's objective, discretized with respect to the state space according to a valid discretization scheme.
   - `funcODE`: Dynamic equations (ODEs)
   - `funcSTAT`: Static equilibrium conditions
   - `Obj`: Objective function (welfare)
   - `mstate`: Nr of backward looking equations
   - `mfor`: Nr of forward looking equations
 You may also set up a new case following the same strucure as the two example cases.


2. Update variable deinitions and classifications:
   - `statevars`: Backward-looking state variables
   - `forwardvars`: Forward-looking jump variables
   - `contvars`: Static control variables

3. Provide an appropriate file to solve for the steady state, and make sure you call it properly in `buildfile.m` (in the current example, `HANK_SSfile.m` or `RANK_SSfile.m`)

4. Re-run `buildfile.m` to regenerate Dynare files

5. Edit the master mod file to include the generated sub-files (in the current example, `HANK.mod` or `RANK.mod`)

### Important things to take into account

- Only one time derivative per equation (must appear linearly)
- Do not mix state and forward variables in the same equation
- Variable names cannot start with `v_` or `dot_`
- Cannot use `pi` as a variable name



## Technical Details

### Discretization

Continuous-time equations are discretized using forward/backward differences:
- State variables: `dot_x ≈ (x - x(-1))/Deltat`
- Forward variables: `dot_x ≈ (x(+1) - x)/Deltat`

### Timing Conventions

- State variables indexed with `(@{Is})`: predetermined
- Forward variables indexed with `(@{If})`: jump variables
- Control variables: contemporaneous (no timing subscript)

## Citation

If you use this code in your research, please cite:

```bibtex
@article{gonzalez2026firm,
  title={Firm Heterogeneity, Capital Misallocation and Optimal Monetary Policy},
  author={Gonz{\'a}lez, Beatriz and Nu{\~n}o, Galo and Thaler, Dominik},
  year={2026}
}
```

## References

Related papers and methods:
- Achdou et al. (2022): "Income and Wealth Distribution in Macroeconomics: A Continuous-Time Approach"
- Kaplan, Moll & Violante (2018): "Monetary Policy According to HANK"
- Auclert et al. (2021): "Using the Sequence-Space Jacobian to Solve and Estimate Heterogeneous-Agent Models"

## License

This project is licensed under the MIT License - see below for details:

```
MIT License

Copyright (c) 2025 Beatriz González, Galo Nuño, Dominik Thaler

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
```


