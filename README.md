# RNAligner - RNA alignement benchmarking tool

RNAligner is a benchmarking tool for state of the art
and legacy RNA alignement prediction algorithms.

## Current algorithms to benchmark and compare:

 - Nussinov
 - ViennaRNA

## Getting an alignement for your first RNA sequence

### 1. Import the necessary structs and functions.

The two structs below allow us to build an RnaSequence structure. This struct can be passed in the Score struct to specify vairous folding algorithms and predict various rna structres.

```rust
use rnaligner::io::RnaSequence;
use rnaligner::compare::Score;
```

### 2. Specify the data that will be used.
 - id: the unique identifier of the sequence. 
    - Note that the id is not really useful in this demonstration, since we only have one sequence to deal with. The id will become more useful when dealing with multiple RNA sequences.

 - seq: the RNA sequence used to predict the folded 2d structure

 - exp_fold: the experimental 2d structure used to evaluate the prediction made by our algorithm.


```rust
let id = "123";
let seq = "AAAUAUGAAGCGAUUUAUUGCAAUUAGUUUCGACCUAAUCUUAGGUGAAAUUCACCCAUAUUUUCCA";
let exp_fold = "(((((((..((((....)))).(((((.......)))))....((((.....)))))))))))....";
```

### 3. Create an RnaSequence struct

```rust
let rna_seq = RnaSequence::new(id, exp_fold, seq); 
```

### 4. Now let's calculate the scores for each algorithm

Here we are running both algorithms to see how they stack up against each other:

```rust
// Calculate the score and clone
let nussinov_score = Score::new(rna_seq.clone(), "nussinov");
let vienna_score = Score::new(rna_seq, "vienna");

// Use the repr() function to print them in the terminal
let _ = nussinov_score.expect("error nussinov").repr();
let _ = vienna_score.expect("error vienna").repr();
```

The scores will tell us how accurate each algorithm is compared to the experimental structure.

You should see something that looks like this:
```
========== Example for id: 123 =============
Algorithm used: nussinov
Experimental result: (((((((..((((....)))).(((((.......)))))....((((.....)))))))))))....
Algorithmic result:  .((((((..(.(((((...).)))(((((((.((((.(...))))))))))).)))))))))(...)
Matches:              ||||||||| ||   |  |    |||    |       ||            ||||||||| ||| 
Match score:         46.27%
Algorithm used: vienna
Experimental result: (((((((..((((....)))).(((((.......)))))....((((.....)))))))))))....
Algorithmic result:  (((((((..((((....))))....((((((.(((((....)))))))))))....)))))))....
Matches:             ||||||||||||||||||||||   ||    |       ||               |||||||||||
Match score:         56.72%

```

### Public associated methods on Score

```rust

pub struct Score
pub fn new(seq: RnaSequence, algo: &str) -> Result<Self, Box<dyn std::error::Error>>
pub fn get_id(&self) -> Result<&str, Box<dyn std::error::Error>>
pub fn get_seq(&self) -> Result<&str, Box<dyn std::error::Error>>
pub fn get_algo(&self) -> Result<String, Box<dyn std::error::Error>>
pub fn get_score(&self) -> Result<String, Box<dyn std::error::Error>>
pub fn get_exp_fold(&self) -> Result<f64, Box<dyn std::error::Error>>
pub fn get_fold(&self) -> Result<f64, Box<dyn std::error::Error>>
pub fn repr(&self) -> Result<(), Box<dyn std::error::Error>>
```

## Running benchmarks on multiple sequences

If you wanna test a bunch of RNA sequences at once (which is way more useful), here's how:

```rust
use rnaligner::io::parse_fasta;
use rnaligner::benchmark::Benchmark;

// Load sequences from a file (395 is the max my machine can handle with Nussinov)
let seq_list: Vec<RnaSequence> = parse_fasta("data/trna_unmodified_dot_bracket.txt", 100);

// Create and run the benchmark
let bench = Benchmark::new(seq_list);
bench.repr();
```

You should see something that looks like this:
```
=========== Benchmark matching scores for VIENNARNA and NUSSINOV algo ==========
Average Match Score for nussinov:  50.26%
Average Match Score for ViennaRNA: 73.04%
Total samples number for each: 100

----------- Min / Max Scores -----------
Nussinov:  min=30% | max=77% | std=11%
ViennaRNA: min=31% | max=100% | std=17%

--- Score Distribution (vienna) ---
0-20%   |   0 | 
20-40%  |   3 | ███
40-60%  |  22 | ██████████████████████
60-80%  |  38 | ██████████████████████████████████████
80-100% |  37 | █████████████████████████████████████

--- Score Distribution (nussinov) ---
0-20%   |   0 | 
20-40%  |  18 | ██████████████████
40-60%  |  61 | █████████████████████████████████████████████████████████████
60-80%  |  21 | █████████████████████
80-100% |   0 | 

Quantity: 100
----- Top five Nussinov -----
1. 76.6% tdbR00000586
2. 76.4% tdbR00000403
3. 73.1% tdbR00000589
4. 70.8% tdbR00000167
5. 70.7% tdbR00000299

----- Worst five Nussinov -----
1. 30.3% tdbR00000093
2. 30.3% tdbR00000097
3. 31.1% tdbR00000553
4. 34.1% tdbR00000410
5. 34.5% tdbR00000592

----- Top five ViennaRNA -----
1. 100.0% tdbR00000089
2. 100.0% tdbR00000493
3. 100.0% tdbR00000053
4. 100.0% tdbR00000080
5. 100.0% tdbR00000015

----- Worst five ViennaRNA -----
1. 31.1% tdbR00000553
2. 37.9% tdbR00000592
3. 38.4% tdbR00000029
4. 41.7% tdbR00000409
5. 44.4% tdbR00000398

```

### Public associated methods on Benchmark

```rust

pub struct Benchmark
pub fn new(seq_list: Vec<RnaSequence>) -> Self
pub fn score_distribution(&self, algo: &str) -> Result<(), Box<dyn std::error::Error>>
pub fn min_max(&self, algo: &str)
pub fn top_five(&self, algo: &str)
pub fn bottom_five(&self, algo: &str)
pub fn repr(&self)
```

## What's going on under the hood?

The tool currently implements two folding algorithms:

1. **Nussinov Algorithm**: This is like the grandpa of RNA folding - simple but gets the job done. It's not super accurate but helps understand the basics.

2. **ViennaRNA**: This is the fancy one. It's wrapped from the C implementation and uses way more sophisticated rules for prediction.

## Building the project on macos

The ViennaRNA package is available in the brewsci/bio tap, which needs to be added to your Homebrew installation.
```bash
brew tap brewsci/bio
```

Once the tap is added, you can install the ViennaRNA Package using the brew install command:
```bash
brew install viennarna
```

Then, run:
```bash
cargo build
cargo run
```

The C code for ViennaRNA will be compiled automatically thanks to the build.rs script.



