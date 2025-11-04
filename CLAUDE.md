# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Summary

The goal of this project is to understand cynipini-plant interactions through dual approaches: analyzing cynipini larvae RNA-seq data (insect side) and plant transcriptomic responses during gall formation (plant side). 

Claude is here to act as a bioinformatic code assistant. It should avoid suggesting biological interpretations unless they indicate clear data quality problems that require technical solutions. It should refrain from judgment about whether a task is adequately completed. The primary goal of our work together is NOT to reach desired output files, but to create a robustly repeatable, documented, and comprehensible workflow so that we can later describe and repeat our methods in collaborations and publications. This is why it is so important not to create hotfixes for code: always make modifications to planned scripts, and avoid proliferation of single-use scripts. Rely on the plan, test, symlink system described below, and document changes to it! Keep scripts with their outputs (ie, create the output folder first, then write the script in that folder, run it from there). 

The species studied are *Andricus quercusfoliatus* (AF), *Callirhytis quercusbatatoides* (CB), *Disholcaspis cinerosa* (DC), *Druon quercuslanigerum* (DL), and *Neuroterus howertoni* (NH). Additional data from *Belonocnema kinseyi* (BK, also labeled as btreatae in some files) are also available from previous work. BK files include both sexual generation (S) and agamic generation (A) larvae. There are also other BK life stages available for comparison. BIP is Biorhiza pallida and DQP is Dryocosmus quercuspalustris. 

## Universal Workflow Best Practices

**IMPORTANT**: Claude Code always starts in the `/carc/scratch/projects/emartins/2016456/adam` directory. The `docs` directory is at `/carc/scratch/projects/emartins/2016456/adam/docs`. Always check your current working directory with `pwd` before looking for relative paths.

Before we begin work today, review and summarize to the user that you understand the following best practices for managing our workflow.

**First Task for New Analyses**: If starting a new analysis pipeline, create a WORKFLOW_BLUEPRINT.md in the project directory that outlines:
- Expected input files and their locations
- Expected output files that will be created
- The sequence of scripts needed
- Key parameters and design decisions

This blueprint serves as the "contract" that all scripts must follow.

**CRITICAL: Always Inspect Input Data Before Analysis**: Before running any analysis pipeline or long-running job:
- Use `head` to examine the first few lines of input files
- Verify the file format (FASTA protein vs nucleotide, TSV column headers, etc.)
- Check file sizes with `ls -lh` to ensure they're reasonable
- For FASTA files, confirm whether sequences are nucleotide (ATGC) or protein (amino acids)
- For tabular data, verify column count and format match expectations
- **Never assume file format from the filename alone** - always verify the actual content

This simple check can prevent hours of wasted computation on incorrectly formatted inputs.

**Project-Specific Status Checking:**
- **Node awareness**: Always check which node you're on with `hostname` and check job queues on both nodes:
  - Local node: `squeue -u $USER`
  - Other node: `ssh hopper 'squeue -u akranz8174'` (if on easley) or `ssh easley 'squeue -u akranz8174'` (if on hopper)
  - Jobs can be submitted cross-node and will run on available resources
  - **NOTE**: Jobs named "claude_dev" are salloc sessions for Claude Code interactive work and can be safely ignored
  - **IMPORTANT**: Do NOT use `head` or `tail` on squeue outputs - they don't work properly. Just use `squeue -u $USER` directly
- For insect work: Check `docs/current_status_insect.md` for job IDs, goals, names, and output locations
- For plant work: Check `docs/current_status_plant.md` for job IDs, goals, names, and output locations
- Then use job IDs to query status with squeue and bash commands that check outputs and log/err files from completed jobs to bring the user up to date on the progress, failures, and successes of each job. 

After conferring with the user, update the appropriate project-specific current_status file accordingly.

## 🗄️ NAS BACKUP & SCRATCH MANAGEMENT

**CRITICAL INODE AWARENESS**: CARC scratch has a strict 1,150,000 file (inode) limit SHARED across the entire project. This is MORE important than storage limits. Trinity assemblies create 100k-170k temp files each. OrthoFinder can create 400k+ files. You MUST manage file counts proactively.

**Complete details**: See `/carc/scratch/projects/emartins/2016456/adam/CARC_SCRATCH_MANAGEMENT.md`

### Quick Reference: NAS Operations

**Check what's backed up on NAS:**
```bash
sshpass -p "WK8y8kcg" sftp akranz8174@nasonia.unm.edu << 'EOF'
cd /NetBackup/carc/scratch/projects/emartins/2016456/adam
ls -l
bye
EOF
```

**Upload directory to NAS (from login node only, NOT SLURM):**
```bash
nohup sshpass -p "WK8y8kcg" sftp akranz8174@nasonia.unm.edu << 'EOF' > backup.log 2>&1 &
cd /NetBackup/carc/scratch/projects/emartins/2016456/adam/destination_directory
put -r source_directory
bye
EOF
```

**Download from NAS to scratch:**
```bash
sshpass -p "WK8y8kcg" sftp akranz8174@nasonia.unm.edu << 'EOF'
cd /NetBackup/carc/scratch/projects/emartins/2016456/adam/source_directory
get -r directory_name
bye
EOF
```

### Standard Workflow for Large Jobs

1. **Before starting**: Check current inode usage and estimate needs
   ```bash
   find /carc/scratch/projects/emartins/2016456/adam -type f 2>/dev/null | wc -l
   ```

2. **Bring only required inputs to scratch** (from NAS if needed)

3. **Run job**

4. **Immediately after completion**:
   - Delete Trinity `read_partitions/` directories (frees 100k-170k inodes each!)
   - Back up final outputs to NAS
   - Verify backup succeeded
   - Delete from scratch

5. **Keep scratch clean**: Only active work belongs on scratch

### Critical Rules

- **ONLY run NAS transfers from login nodes** (easley/hopper), NEVER via SLURM
- **Delete Trinity temp files immediately** after assembly completes
- **Verify backups before deleting** from scratch
- **Coordinate with project members** (daniel.paulo, soham) before large jobs

## 📐 WORKFLOW BLUEPRINT SYSTEM

**CRITICAL**: Before starting any analysis, create a WORKFLOW_BLUEPRINT.md that defines:
- **Input files**: Exact paths and expected formats
- **Output files**: Exact paths that downstream scripts will expect
- **Processing steps**: Numbered scripts with clear purposes
- **Dependencies**: Use Beads (`bd`) to track step dependencies for complex workflows

**For complex multi-step workflows**: Initialize Beads in the project directory:
```bash
cd project_directory/
bd init
bd create "Step 1: Parse input data" --tag=step1
bd create "Step 2: Run analysis" --tag=step2
bd dep add step2 step1  # step2 depends on step1
bd list
```

**Scripts live WITH their outputs** - no separate scripts/ folder. Each processing step gets its own directory with numbered scripts.

## 🔄 VERSIONED OUTPUT WITH SYMLINKS

**All scripts must follow this pattern**:

1. **Versioned outputs**: Include timestamp in output filename
   ```python
   timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
   output_file = f"comparison_at_matrix_{timestamp}.tsv"
   ```

2. **Update symlink to current**: Point to latest successful output
   ```python
   if success:
       symlink = "comparison_at_matrix_CURRENT.tsv"
       if os.path.exists(symlink):
           os.unlink(symlink)
       os.symlink(output_file, symlink)
   ```

3. **Downstream scripts read from CURRENT**:
   ```python
   matrix = pd.read_csv("comparison_at_matrix_CURRENT.tsv")
   ```

## 📝 PROVENANCE TRACKING

**Every script must output provenance**:
```python
provenance = {
    'script': os.path.abspath(__file__),
    'timestamp': datetime.now().isoformat(),
    'inputs': {
        'cohort_file': os.path.abspath(cohort_file),
        'de_results': [os.path.abspath(f) for f in de_files]
    },
    'outputs': {
        'matrix': os.path.abspath(output_file),
        'symlink': os.path.abspath(symlink)
    },
    'parameters': {
        'min_genes': min_genes,
        'aggregation_method': 'weighted_mean'
    },
    'git_hash': subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode().strip()
}
with open(f"provenance_{timestamp}.json", 'w') as f:
    json.dump(provenance, f, indent=2)
```

## 🧩 MODULAR SCRIPT DESIGN

**For time-consuming analyses, design for SLURM arrays**:

```python
# Instead of:
for study in all_studies:
    process_study(study)

# Design as:
if __name__ == "__main__":
    study_id = sys.argv[1]  # Or use SLURM_ARRAY_TASK_ID
    process_study(study_id)
```

```bash
#SBATCH --array=0-21  # For 22 studies
studies=(study_01 study_02 ... study_22)
python process_study.py ${studies[$SLURM_ARRAY_TASK_ID]}
```

## 📊 INODE-CONSCIOUS WORKFLOW DESIGN

**CRITICAL**: Design workflows to minimize file proliferation. The 1.15M inode limit is MORE constraining than storage space.

### High-Inode Operations to Avoid

**❌ BAD: Creating many intermediate files**
```python
# Creates 10,000+ small files
for gene_id in gene_list:
    with open(f"intermediate/{gene_id}_result.txt", 'w') as f:
        f.write(process_gene(gene_id))
```

**✅ GOOD: Aggregate into single files**
```python
# Creates 1 file with all results
results = {gene_id: process_gene(gene_id) for gene_id in gene_list}
with open("all_results.json", 'w') as f:
    json.dump(results, f)

# Or use pandas for tabular data
pd.DataFrame(results).to_csv("all_results.tsv", sep='\t')
```

### Tools Known to Create Many Files

**Trinity**: 100k-170k temp files in `read_partitions/`
- **Solution**: Delete immediately after assembly: `rm -rf trinity_out/read_partitions/`
- Keep only: `Trinity.fasta`, `Trinity.fasta.gene_trans_map`

**OrthoFinder**: 400k+ files for large proteome sets
- **Solution**: Back up to NAS immediately, delete from scratch
- Only bring back summary files for analysis

**Per-gene/per-transcript outputs**: AVOID whenever possible
- **Instead**: Use multi-sequence FASTA files or tabular formats
- **Example**: Don't create 50k individual protein FASTA files; use one multi-FASTA

**SignalP/DeepLoc batch outputs**: Can create 100k+ individual result files
- **Solution**: Aggregate results into single TSV immediately
- Delete individual output files after aggregation

### Design Principles

1. **Consolidate early**: Merge per-sample/per-gene results into summary tables ASAP
2. **Use archives for many small files**: If you must keep many small files, tar them:
   ```bash
   tar -czf results.tar.gz results_directory/
   rm -rf results_directory/
   ```
3. **Clean up intermediate files in pipelines**: Add cleanup steps to SLURM scripts
4. **Estimate inode cost before running**:
   ```bash
   # Check current count
   find /carc/scratch/projects/emartins/2016456/adam -type f | wc -l
   # Estimate job impact (e.g., 5 Trinity assemblies = ~650k inodes)
   ```

## 🚫 WORKFLOW ANTI-PATTERNS TO AVOID

1. **NEVER create patch scripts**: If study_22 is missing, fix the main script to include it, don't write `add_study22.py`
2. **NEVER use ambiguous names**: No `_fixed`, `_corrected`, `_final`, `_new`
3. **NEVER modify outputs in place**: Always create new versioned output + update symlink
4. **NEVER break downstream paths**: Scripts should read from `_CURRENT` symlinks
5. **NEVER proliferate small files**: Avoid creating thousands of per-gene/per-sample files; use aggregated formats (multi-FASTA, TSV, JSON)

## 🎯 BLUEPRINT-FIRST DEVELOPMENT

**CRITICAL**: WORKFLOW_BLUEPRINT.md is the single source of truth. Always check and update it before implementing.

**Key Principles**:
- ONE blueprint per analysis
- ONE location for each output type
- ONE current version of each script
- Fix the authoritative script, not the output
- Update blueprint when adding functionality

**When fixing missing data** (e.g., study_22):
1. Check blueprint for where data should appear
2. Fix the script that creates that output
3. Run script → new versioned output → symlink updates automatically  

## Environment Management

**Micromamba**: Always use micromamba, not conda, for new environments. Available environments listed in `docs/micromamba_library.md`.

**CRITICAL PYTHON ENVIRONMENT REQUIREMENT**:
- **ALWAYS activate busco_env before running any Python scripts**: Use `micromamba activate busco_env` before any Python commands
- This environment contains essential packages like pandas, numpy, matplotlib, biopython, etc.
- Never run Python scripts without first activating this environment - they will fail with "ModuleNotFoundError"

**Go Language**: Installed at `/carc/scratch/projects/emartins/2016456/tools/go`
- Go 1.22.0 is available for building Go-based tools
- Environment variables set in `~/.bashrc`:
  - `PATH` includes `/carc/scratch/projects/emartins/2016456/tools/go/bin`
  - `GOPATH=/carc/scratch/projects/emartins/2016456/go`
  - `$GOPATH/bin` added to PATH for Go-installed tools

**Beads Task Manager** (`bd`): Lightweight issue tracker with dependency support
- Installed at `$GOPATH/bin/bd`
- Use for tracking complex multi-step workflows
- **Organization**: One Beads DB per major analysis directory (e.g., `output2/meta_clustering/.beads/`)
  - Each project gets its own `bd init` in its root directory
  - Keeps task tracking co-located with the work, matching the "scripts with outputs" philosophy
  - Navigate to project directory before running `bd` commands
- See `bd quickstart` for usage or https://github.com/steveyegge/beads

**CRITICAL SEQUENCE LABELING RULE**: When extracting sequences from any individual Trinity assembly, ALWAYS append a species/sample label indicating its origin (e.g., AFL32_TRINITY_DN123, DQP_TRINITY_DN456, BIP_TRINITY_DN789). Never create sequence datasets without proper origin labeling as this makes downstream phylogenetic analyses impossible to interpret.

Assume that the project-specific data_locations files (data_locations_insect.md or data_locations_plant.md) provide locations to important files not identified in recent context; don't use directory-wide searches by keyword to find existing files. If you can't find it, ask the user rather than searching. 

When problems occur, check docs/programming_issues.md for previous solutions and troubleshooting.

As needed, look up important project information in the following files: 

**Universal files:**
- docs/micromamba_library.md contains a list of existing micromamba envs
- docs/programming_issues.md lists solutions to previous frequently encountered coding issues with key packages and software

**Project-specific files:**
- docs/current_status_insect.md and docs/current_status_plant.md contain current job status and progress
- docs/data_locations_insect.md and docs/data_locations_plant.md contain important pathways for frequently used, validated datasets

SLURM job info:
Start scripts with #!/bin/bash
Max memory: 90GB
CPUs max: 64
# nodes: 2
Processors per node: 32
Max time limit: 2 days
Partition: general (default)

**CRITICAL: SLURM jobs start in the directory where sbatch is called, NOT where the script lives!**
- **ALWAYS add `cd /full/path/to/project/directory` at the start of every SLURM script**
- This prevents scripts from creating outputs in the home directory when submitted via ssh
- Example:
```bash
#!/bin/bash
#SBATCH directives...

cd /carc/scratch/projects/emartins/2016456/adam/output2/my_project
# Now run commands...
```

IMPORTANT: Use SLURM's printf-style placeholders for output/error files, not shell variables. `$SLURM_JOB_ID` isn't set when `#SBATCH` lines are parsed.

For single jobs:
```
#SBATCH --output=my_job_%j.out
#SBATCH --error=my_job_%j.err
```

For array jobs:
```
#SBATCH --array=0-9
#SBATCH --output=my_job_%A_%a.out
#SBATCH --error=my_job_%A_%a.err
```

Useful tokens: `%j`=jobid, `%A`=array master id, `%a`=array index, `%x`=job name, `%u`=user, `%N`=node

**NEVER USE THE `find` BASH COMMAND** - Use list, LS, Grep, or other tools instead. Only use `find` if it's the ONLY way to accomplish a specific task and no other tool will work.

**NEVER USE `head`, `tail`, or `grep` ON `squeue` OUTPUT** - The squeue command output doesn't work properly with pipes. Always use `squeue -u $USER` by itself without any additional processing.

**ALWAYS USE FULL PATHS WHEN SUBMITTING JOBS ACROSS NODES** - When using `ssh easley` or `ssh hopper` to submit jobs, always use the full path to the script (e.g., `ssh easley "sbatch /carc/scratch/projects/emartins/2016456/adam/script.slurm"`), not relative paths.

**NEVER USE GLOBAL SEARCH (Grep/Glob without specific paths) - IT CRASHES THE TERMINAL!** - Always start by using `ls` to examine the directory structure first. Make intelligent choices about which specific directories to search within based on what you see, or ask the user for guidance. Never run broad searches across the entire project directory tree - this will crash your terminal session and lose all progress.

## 🔧 AVOIDING GLOB ERRORS IN BASH
When using wildcards with pipes or redirects, bash fails to expand globs properly, causing "cannot access 'glob'" errors.

**ROOT CAUSE**: The Bash tool doesn't properly handle glob expansion when combined with pipes or redirects (2>/dev/null).

**THE PROBLEM**:
```bash
# FAILS - produces "ls: cannot access 'glob': No such file or directory"
ls *.txt | head -5
ls *.txt 2>/dev/null | head -5
find /path -name "*.log" 2>/dev/null | head

# Also FAILS with brace expansion:
ls /path/salmon*.{err,out,log} 2>/dev/null
```

**SOLUTION 1 - Use bash -c (MOST RELIABLE)**:
```bash
# Wrap the entire command in bash -c to ensure proper glob expansion
bash -c "ls *.txt | head -5"
bash -c "ls /path/salmon*.{err,out,log} 2>/dev/null | head -5"
```

**SOLUTION 2 - Use find without stderr redirect**:
```bash
# For finding files, use find with proper parentheses
find /path -maxdepth 2 \( -name "*.out" -o -name "*.err" \) | head -10

# NOT: find /path -name "*.out" 2>/dev/null | head
```

**SOLUTION 3 - Use nullglob for arrays**:
```bash
shopt -s nullglob
files=(*.txt)
if [ ${#files[@]} -gt 0 ]; then
    printf '%s\n' "${files[@]}" | head -5
fi
```

**SOLUTION 4 - Avoid globs with xargs**:
```bash
# Use find with xargs instead of glob expansion
find . -maxdepth 1 -name "*.txt" -print0 | xargs -0 ls -la | head -5
```

**Why this happens**: The Bash tool's shell environment doesn't properly expand globs before parsing pipes and redirects, treating the glob pattern as a literal string.

## Work Scoping and Planning

**For complex multi-step tasks**: Use Beads (`bd`) to track phases and dependencies:
```bash
bd create "Phase 1: Sequence prep" --tag=phase1
bd create "Phase 2: Tree building" --tag=phase2
bd create "Phase 3: Visualization" --tag=phase3
bd dep add phase2 phase1
bd dep add phase3 phase2
```

**NEVER use placeholder content** in scripts or outputs. Complete each phase fully before moving to the next.

## File Organization and Documentation Strategy

**CRITICAL FILE ORGANIZATION RULES**: Follow the versioned output system with workflow blueprints.

### Project-Specific File Organization

**CRITICAL: NO SEPARATE SCRIPTS FOLDER!**

Scripts must live WITH their outputs in the same directory. This prevents script archaeology and ensures clear relationships between code and results.

**Key Principles:**
- **Scripts live in the SAME directory as their outputs** (never in scripts/)
- Each processing step gets its own directory
- Scripts are numbered to match workflow steps
- Versioned outputs accumulate (don't delete old versions)
- CURRENT symlinks always point to latest successful output
- Every output has a corresponding provenance JSON

### Script Versioning

**Version Scripts Chronologically**:
```bash
# Good naming:
aggregate_to_at_v1.py          # First attempt
aggregate_to_at_v2.py          # Fixed bug
aggregate_to_at_v3.py          # Added study_22

# Also acceptable:
aggregate_to_at_20250812.py    # Date-based
aggregate_to_at_3f2a1b.py      # Git commit hash
```

**NEVER use**:
- `script_fixed.py`, `script_corrected.py`, `script_final.py`
- These names provide no chronological information

### Documentation Workflow

**Current Status Files** (`docs/current_status_[side].md`):
- Document ongoing work and active job status
- Record outcomes (what worked vs failed)

**For complex workflows**: Use Beads to track status, document outcomes in current_status files.

**Transfer to Historical Record**: ONLY after user confirms success. Never preemptively mark as "complete".

### File Cleanup Protocol

**After Each Session:**
- Move .out/.err files to appropriate `logs/` folder
- Organize scripts into project-specific folders
- Move temporary/debug scripts to `temp_scripts/`
- Clean up failed analysis attempts

**User Approval Required Before:**
- Marking any analysis as "complete" in documentation
- Moving workflows from current_status to historical records
- Deleting any files (always ask first)

## 🧪 TESTING AND VALIDATION FRAMEWORK

**CRITICAL**: All new pipelines must include validation from the start. Test-driven development prevents the "fix-and-fork" cycle.

### Core Testing Principles

1. **Data Contracts**: Define explicit schemas for all inputs/outputs
   ```python
   import pandera as pa
   
   HitsSchema = pa.DataFrameSchema({
       "study": pa.Column(str, nullable=False),
       "qseqid": pa.Column(str, nullable=False), 
       "score": pa.Column(float, checks=pa.Check.greater_than(0))
   })
   
   @pa.check_io(df=HitsSchema)
   def process_hits(df: pd.DataFrame) -> pd.DataFrame:
       # Function will validate input/output automatically
   ```

2. **Pure Functions**: Separate I/O from logic
   ```python
   # BAD: Mixed concerns
   def process_data(filename):
       df = pd.read_csv(filename)  # I/O mixed with logic
       return df.groupby('study').mean()
   
   # GOOD: Pure function
   def process_data(df: pd.DataFrame) -> pd.DataFrame:
       return df.groupby('study').mean()  # Pure logic, testable
   ```

3. **Property-Based Testing**: Test invariants, not specific values
   ```python
   # Test: "All input natives get assigned or explicitly dropped"
   def test_coverage_invariant(input_df, output_df, drop_log):
       assert len(output_df) + len(drop_log) == input_df['native_id'].nunique()
   ```

4. **Deterministic Sampling**: For test data, use stable hashing
   ```python
   def sample_for_testing(df: pd.DataFrame, fraction: float = 0.01):
       df['hash'] = df['native_id'].apply(lambda x: hash(x) % 10000)
       return df[df['hash'] < fraction * 10000]
   ```

### Required Validation Outputs

Every pipeline must produce:

1. **Coverage Report**: What percentage of inputs produced outputs?
2. **Drop Log**: Why did specific records fail? (with counts by reason)
3. **Sanity Checks**: Basic assertions (e.g., no duplicate IDs, all studies present)

### Test Structure for New Analyses

```
my_analysis/
├── WORKFLOW_BLUEPRINT.md
├── src/
│   ├── 01_process_data.py      # Main pipeline script
│   └── core/                   # Pure functions (no I/O)
│       └── assignment.py
├── tests/
│   ├── test_data/             # Small test datasets
│   ├── test_core.py           # Unit tests for pure functions
│   ├── test_integration.py    # Full pipeline tests
│   └── test_validation.py     # Output validation tests
└── validate.py                # Run all validations
```

### Validation Script Template

Every pipeline should include a validation script:

```python
#!/usr/bin/env python3
"""validate.py - Run after pipeline completes"""

def validate_pipeline_outputs(output_dir):
    results = {}
    
    # Check all expected files exist
    expected_files = ["assignments.json", "matrix.tsv", "drop_log.json"]
    for f in expected_files:
        results[f] = os.path.exists(output_dir / f)
    
    # Check coverage
    coverage = calculate_coverage(input_file, output_file)
    results['coverage'] = coverage
    assert coverage > 0.95, f"Coverage {coverage} below 95% threshold"
    
    # Check no studies lost
    studies_in = get_unique_studies(input_file)
    studies_out = get_unique_studies(output_file)
    results['missing_studies'] = studies_in - studies_out
    assert not results['missing_studies'], f"Lost studies: {results['missing_studies']}"
    
    return results
```

### When to Write Tests

1. **Before starting**: Write validation criteria
2. **When debugging**: Add test that reproduces the bug
3. **After fixing**: Ensure all tests pass

This testing framework is mandatory for all new analyses to ensure reliability and reproducibility.

## 🔧 GIT VERSION CONTROL FOR SOFTWARE PIPELINES

**CRITICAL NFS ISSUE**: Git databases (.git directories) CANNOT be stored on /carc/scratch NFS filesystem - they become corrupted. Use separated git directory approach.

### Setting Up Git for Pipeline Development

**For long-term software pipelines** (multi-script workflows that will be published/shared):

1. **Store git database in home directory** (stable filesystem):
   ```bash
   cd /path/to/pipeline/on/scratch
   git init
   mv .git /users/akranz8174/.git-repos/pipeline_name
   echo "gitdir: /users/akranz8174/.git-repos/pipeline_name" > .git
   ```

2. **Configure git**:
   ```bash
   git config user.name "Adam Kranz"
   git config user.email "adamjameskranz@gmail.com"
   git branch -m main
   ```

3. **Test with a commit**:
   ```bash
   git add .gitignore
   git commit -m "Initial commit"
   ```

### What to Commit

**DO commit**:
- Python/R/Shell scripts (numbered pipeline scripts)
- WORKFLOW_BLUEPRINT.md and documentation
- Configuration files (config.py, etc.)
- Small input files (<100KB) like gene lists, metadata tables
- .gitignore file

**DO NOT commit**:
- Output files (.tsv, .fasta, .json results)
- Log files (.out, .err, .log)
- Data directories (raw reads, assemblies, databases)
- Test outputs and temporary files
- Python cache (__pycache__/)
- Large reference files (>100KB)

### .gitignore Template for Bioinformatics Pipelines

```gitignore
# Outputs and data
outputs/
logs/
data/
inputs/
test_outputs/
*.out
*.err
*.log

# Results files
*.fasta
*.fastq
*.bam
*.sam

# Python cache
__pycache__/
*.pyc

# Temporary and test files
temp_*
test_*

# Beads tracker (project-specific)
.beads/

# Keep templates but ignore specific data
*_CURRENT.tsv
*_[0-9][0-9][0-9][0-9][0-9][0-9]*.tsv
```

### Commit Message Format

**For pipeline development**, use descriptive commit messages:

```bash
# Good commit messages:
git commit -m "Add frame-aware clustering to Phase 4 target detection"
git commit -m "Fix: Gene grouping now correctly combines exons"
git commit -m "Add adaptive windowing strategy (0-200kb) to Phase 6"

# Document what's missing:
git commit -m "Phase 6 fixes complete
- Frame-aware clustering
- Adaptive windowing
- Gene grouping

MISSING: Tandem LOC deduplication for batch processing"
```

### Connecting to GitHub

1. **SSH key already configured** at `~/.ssh/id_ed25519.pub`
2. **Create repo on GitHub**: https://github.com/new
3. **Add remote and push**:
   ```bash
   git remote add origin git@github.com:Megachile/repo_name.git
   git push -u origin main
   ```

### When to Use Git

**Use Git when**:
- Building a multi-script pipeline for publication
- Multiple people will work on the code
- Pipeline will be maintained long-term
- Code changes frequently and you need version history

**Don't use Git for**:
- Single-use analysis scripts
- Quick data explorations
- Scripts that only produce outputs (use versioned outputs + symlinks instead)
- Anything where the workflow blueprint + provenance JSON is sufficient

### Workflow: Making Changes

```bash
# 1. Make changes to scripts
vim scripts/04_blast_targets.py

# 2. Test the changes
sbatch run_phase4.slurm

# 3. Commit when working
git add scripts/04_blast_targets.py
git commit -m "Fix: Add missing frame validation in target clustering"

# 4. Push to GitHub
git push
```

### Force Push (When Safe)

**Force push is SAFE when**:
- You're the only developer
- Repo is brand new (first few commits)
- You need to clean up history

```bash
git push --force origin main
```

**Force push is DANGEROUS when**:
- Others are using the repo
- Repo is shared/published
- You're unsure what you're overwriting