# data-workflows-tutorial

[Snakemake](https://snakemake.readthedocs.io/en/stable/) data workflow tutorial repository.

Snakemake is a workflow manager that is strongly coupled with Python.
If you are comfortable with scripting in Python, it shouldn't take long to become familiar with Snakemake!

Workflow managers allow you to define a distinct set of data processing steps without worrying about making sure you run them in the correct order and ensuring steps are only re-run if something has changed.
By modularising your thinking on data workflows into these steps, you can keep yourself sane when processing lots of raw data sources and running a simulation/optimisation tool across a variety of scenarios.

## Prepare

> [!TIP]
> If you already have a setup on your device to install conda environments, then feel free to skip ahead to the [github repository cloning section](#clone-the-github-repository).

### Download supporting software

First, you will need to install software on your device:

1. [VSCode](https://code.visualstudio.com/download).
This gives you access to an Interactive Development Environment (IDE) in which to edit code and to interact with your device's terminal.
On Windows, I recommend you follow the _command prompt_ instructions.
2. [Miniforge](https://github.com/conda-forge/miniforge?tab=readme-ov-file#download).
This gives you access to `conda` in your device's terminal, with which you can create isolated Python environments to work in.
There are other variants of access to `conda` (Anaconda, Miniconda, Mambaforge, etc.).
I recommend `Miniforge` as it defaults to the open-source `conda-forge` [channel](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-channels.html), instead of the commercial `defaults` channel.
> [!NOTE]
> On Windows, you should download `Miniforge3-Windows-x86_64` from the releases.
> You may find that on Chrome it says the download is "dangerous" - you can safely bypass this and force the download.

### Set up VSCode

With this software installed, you can then set up VSCode to allow access to `conda` and `git` from the "terminal":

1. **Set the correct terminal**
    1. Open VSCode.
    1. Open the `command palette` (Control-Shift-P, shows a bar at the top of the screen).
    1. Search for `Terminal: Select Default Profile`.
    1. Select `Command Prompt`.
1. **Install the Python Extension**.
Click on the "Extensions" tab on the left-hand side of VSCode (four squares) and search for and install the official Microsoft Python extension.
1. **Set your Python environment**
    1. Open the `command palette` (Control-Shift-P, shows a bar at the top of the screen).
    1. Search for `Python: Select Interpreter`
    1. Select `base` (it should have `miniforge3` in the file path)
1. Open a new terminal window (`Terminal` top-bar tab -> `New Terminal`, it will open at the bottom of the screen).
1. Run `conda install git` in the terminal.

### Clone the GitHub repository

With VSCode set up, you can "clone" (i.e. copy) this repository to your device:

1. In the VSCode terminal, call `git clone https://github.com/brynpickering/data-workflow-tutorial.git <output-directory>\data-workflow-tutorial`.
`<output-directory>` should be a directory on your device where you want to store cloned GitHub repositories (often something like `%USERPROFILE%\Repositories` on Windows or `~/Repositories` on Linux/MacOS).
1. Then you can open that cloned repository (i.e. downloaded folder) in VSCode (`File` top-bar tab -> `Open Folder` -> navigate to `<output-directory>\data-workflow-tutorial`).

### Create the masterclass conda environment

Finally, you can create the isolated Python working environment for this masterclass:

1. In your `<output-directory>\data-workflow-tutorial` VSCode session, open a new terminal.
1. Run the following two commands to create your working environment and then to _activate_ it:
    ```bash
    conda env create -f environment.yaml
    conda activate data-workflow-tutorial
    ```

### Activate the Snakemake Language extension

Click on the "Extensions" tab on the left-hand side of VSCode (four squares) and search for and install the [`Snakemake Language`](https://marketplace.visualstudio.com/items?itemName=Snakemake.snakemake-lang) extension.
This will make it easier to read the workflow files as it will provide highlighting support.

## Understanding the repository structure

Snakemake workflows are composed of a number of key sections / files:

### The Snakefile

This is the powerhouse of a workflow, it collects all other workflow files you define and usually defines an `all` rule, which is the one rule to rule them all.
If you have a small workflow, you may find that you do not need separate workflow files and can instead define all rules inside this file.

### rule files (rules/*.smk)

You will usually have lots of rules that you need to separate out into different files.
You can think of these as "modules" which are used to complete a subset of tasks.

### Isolated environment definitions (envs/*.yaml)

A key component of reproducibility is running code in isolated environments, in which you have strictly pinned dependencies.
These environments will be called upon whenever a rule refers to them with the `conda:` option.
That rule will then be triggered inside that isolated environment.

### Scripts (scripts/*.py)

You will be mostly processing data using Python scripts.
Ideally you have one script linked to one Snakemake `rule`.
Scripts should be small and focussed.
If a script is getting long, it's a sign that you need to split it out into separate scripts with an intermediate dataset that is passed between rules.

### User configuration (config/*.yaml)

You will likely have a large number of resources for and assumptions you apply to your data workflow.
It's useful to define these outside your script and rule files so that they are (a) transparently defined and (b) can be easily updated at a later date if data sources update or your assumptions change.

>[!NOTE]
>`config/default.yaml` is a special case file. It will be loaded automatically when you run your workflow.
> You can define other configuration files, which will be applied _on top of_ the default file (adding to / replacing content in that file).
> You will just need to remember to refer to those files when you run your workflow `snakemake --configfile config/my-config.yaml`

## Create your workflow

To create your workflow, you will need a collection of rule-script-env-io-config combinations.
That is, you need a Snakemake `rule` which will run a `script` inside a conda `env` to process `inputs` into `outputs` with a number of `config` parameters.

For instance, let's say we want to open a CSV file, apply arithmetic to its contents by a configured value and then save it to parquet file then we'll need the following:

### rule

In your Snakefile or a separated `rules/*.smk` file, you need a rule:

```py
rule my_new_rule:
    message: "This is a readable message that will be printed to the terminal when this rule is activated"
    input:
        to_convert=config["resources"]["convert_me"],
    output:
        converted="./results/converted_table.parquet", # use `../resources/<filepath>` if you are writing this in a separate rule file in `rules/`.
    params:
        multiplier=config["parameters"]["multiply_with_me"],
        adder=config["parameters"]["add_with_me"],
    conda:
        "./envs/table_processor.yaml" # use `../envs/<filepath>` if you are writing this in a separate rule file in `rules/`.
    script:
        "./scripts/multiply_data.py"
```

### script

When Snakemake takes control of your process, it will inject the `snakemake` object into your script, which means you can access the contents of the [rule](#rule).
This means you can set up your script to work with the rule data as follows:

```py
# "./scripts/multiply_data.py"
import pandas as pd

input_df = pd.read_csv(snakemake.input.to_convert)
output_df = input_df * snakemake.params.multiplier + snakemake.params.adder
output_df.to_parquet(snakemake.output.converted)
```

We tend to split out the `snakemake` object from the actual data processing so that we could feasibly test our function separately.
We do so by using `if __name__ == "__main__"` so that the snakemake object is only accessed when the script is _executed_ not _imported_.
This means our script becomes:

```py
# "./scripts/multiply_data.py"
import pandas as pd

def run_me(input_df: pd.DataFrame, multiplier: float | int, adder: float | int) -> pd.DataFrame:
    """This is a script that applies arithmetic to a DataFrame.

    Args:
        input_df (pd.DataFrame): DataFrame to which arithmetic will be applied.
        multiplier (float | int): the input values will be multiplied by this value.
        adder (float | int): this value will be added to the input _after_ `multiplier` has been applied.
    Returns:
        pd.DataFrame: table with the same shape as `input_df` but with arithmetic applied.
    """
    output_df = input_df * multiplier + adder
    return output_df

if __name__ == "__main__":
    input_df = pd.read_csv(snakemake.input.to_convert)
    output_df = run_me(input_df, snakemake.params.multiplier, snakemake.params.adder)
    output_df.to_parquet(snakemake.output.converted)
```

### Configuration

We have requested some configuration parameters to run our rule.
We access them from a [YAML](https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html) file - a human-readable Python dictionary.

In this example, our configuration file would look like:

```yaml
# "./config/default.yaml"
resources:
  convert_me: resources/dummy_table.csv
parameters:
  multiply_with_me: 10
  add_with_me: 5.5
```

### Isolated environment

To run this example we need Pandas and the PyArrow engine (to store the data in `.parquet` format).
So we add these to our isolated environment:

```yaml
# "./envs/table_processor.yaml"
name: table_processor
channels:
  - conda-forge  # we ❤️ the `conda-forge` Anaconda channel!
  - nodefaults  # we do not like the `defaults` Anaconda channel.
dependencies:
  # dependencies are strictly pinned so that we know in we'll get the same functionality every time we run in this environment
  - pandas = 2.2.3
  - pyarrow = 19.0.1
```

### Chaining rule-script-env-io-config combinations together

The benefit of data workflow management is that you can chain processes together.

Continuing our example, you could use the output of the rule as the input to another rule:

```py
rule my_second_rule:
    message: "This is a readable message that will be printed to the terminal when this rule is activated"
    input:
        to_process_further="./results/converted_table.parquet", # use `../resources/<filepath>` if you are writing this in a separate rule file in `rules/`.
    output:
        processed="./results/further_processed_table.parquet"
    conda:
        "./envs/table_processor.yaml" # use `../envs/<filepath>` if you are writing this in a separate rule file in `rules/`.
    script:
        "./scripts/further_process_data.py"
```

In this case, we aren't using any parameters (no `params` defined) but we have changed our script reference to a script that will process the table differently.

## Run your workflow

>[!TIP]
>You can check what the workflow would do _without_ actually running it by using the `--dry-run` / `-n` flag.

To run your workflow, you have three options:

### Request a specific output file

To request that your workflow is run just to produce a file of interest, explicitly mention that file:

```sh
snakemake "./results/converted_table.parquet"
```

You can provide multiple files in a list as well:

```sh
snakemake "./results/converted_table_1.parquet" "./results/converted_table_3.parquet"
```

### Request rules be run

>[!WARNING]
>We don't recommend you actually use this, as it will fail as soon as you introduce [wildcards](#wildcards).

From the previous examples, we could request that `my_second_rule` is run.
Once complete, we will expect to see both `converted_table.parquet` and `further_processed_table.parquet` in our results directory.

```sh
snakemake my_second_rule
```

### Run `all`

If you just call `snakemake` it will trigger the `all` rule inside the `Snakefile`.
You can then define your `all` rule to request a number of inputs.
Snakemake will then ensure these files are all generated:

```python
rule all:
    message: "Collecting all required project results"
    input:
        "./results/converted_table_1.parquet",
        "./results/converted_table_2.parquet",
        "./results/converted_table_3.parquet",
```

## Advanced functionality

There is a lot of advanced functionality.
Here, we'll look at a few key things.
If you're having fun then I recommend you spend more time with the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/#).

### Wildcards

[Wildcards](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards) allow you to generalise rules so they can be triggered to create similar types of data.
For example, consider the above example [rule](#rule).
Let's say we want to create lots of files with different multipliers, we can use a wildcard using curly brackets:

```py
rule my_new_rule:
    message: "Apply arithmetic, including a {wildcards.multiplier}x multiplier, to a data table."
    input:
        to_convert=config["resources"]["convert_me"],
    output:
        converted="./results/converted_table_{multiplier}.parquet", # use `../resources/<filepath>` if you are writing this in a separate rule file in `rules/`.
    params:
        adder=config["parameters"]["add_with_me"],
    conda:
        "./envs/table_processor.yaml" # use `../envs/<filepath>` if you are writing this in a separate rule file in `rules/`.
    script:
        "./scripts/multiply_data.py"
```

This wildcard is now available to our script:

```py
# "./scripts/multiply_data.py"
import pandas as pd

def run_me(input_df: pd.DataFrame, multiplier: float | int, adder: float | int) -> pd.DataFrame:
    """This is a script that applies arithmetic to a DataFrame.

    Args:
        input_df (pd.DataFrame): DataFrame to which arithmetic will be applied.
        multiplier (float | int): the input values will be multiplied by this value.
        adder (float | int): this value will be added to the input _after_ `multiplier` has been applied.
    Returns:
        pd.DataFrame: table with the same shape as `input_df` but with arithmetic applied.
    """
    output_df = input_df * multiplier + adder
    return output_df

if __name__ == "__main__":
    input_df = pd.read_csv(snakemake.input.to_convert)
    output_df = run_me(
        input_df,
        float(snakemake.wildcards.multiplier), # <- grabbing the wildcard.
        snakemake.params.adder
    )
    output_df.to_parquet(snakemake.output.converted)
```

Now, we can request files with the wildcard filled in and it will extract the value and pass it to `wildcards.multiplier`.
For instance, we could call:

```sh
snakemake "./results/converted_table_1.parquet" "./results/converted_table_3.parquet"
```

The first file will have multiplied all values in the table by `1`, the second file will have multiplied them by `3`.

We could also use it as an input to another rule:

```py
rule my_second_rule:
    message: "This is a readable message that will be printed to the terminal when this rule is activated"
    input:
        to_process_further_1="./results/converted_table_1.parquet",
        to_process_further_2="./results/converted_table_2.parquet",
        to_process_further_3="./results/converted_table_3.parquet",
    output:
        processed="./results/further_processed_table.parquet"
    conda:
        "./envs/table_processor.yaml" # use `../envs/<filepath>` if you are writing this in a separate rule file in `rules/`.
    script:
        "./scripts/further_process_data.py"
```

This will trigger the `my_new_rule` rule _three times_ to create each of the required inputs to this rule.
Each time, it will pass a different multiplier.

You can use multiple wildcards in a single rule, but every one must be part of the output filename(s).

### Running shell commands

You can do more than call Python scripts in your rules, you can also choose to run a shell command, as you would in the terminal.
For instance, we often want to download data from the internet as a first step in our project:

```py
# Note how there is no `input` as we are grabbing our input from outside the project.
rule download_data:
    message: "Download a table from the UK govt. website"
    params:
        url = config["resources"]["govt_energy_demand_statistics"]
    output:
        downloaded="./resources/msoa_domestic_elec.xlsx"
    conda:
        "./envs/shell.yaml"
    shell:
        "curl -sSLo {output.downloaded} '{params.url}'"
```

In our configuration, our URL would be given as <https://assets.publishing.service.gov.uk/media/6762d1e1cdb5e64b69e3077a/MSOA_domestic_elec_2010-2023.xlsx>.
We have also referenced a different conda environment, as now we need access to the `curl` tool.

### Python code in rules

Snakefiles accept ad-hoc Python code.
This means that you can e.g. define `if` statements or `for` loops inside your rule.

```py
rule my_second_rule:
    message: "This is a readable message that will be printed to the terminal when this rule is activated"
    input:
        *[f"./results/converted_table_{val}.parquet" for val in [1, 2, 3]], # snakemake also has the `expand` helper for this: `expand(./results/converted_table_{val}.parquet", val=[1, 2, 3])`
    output:
        processed="./results/further_processed_table.parquet"
    params:
        val=config["parameters"]["foo"] if config["parameters"]["bar"] == 5 else config["parameters"]["baz"]
    conda:
        "./envs/table_processor.yaml"
    script:
        "./scripts/further_process_data.py"
```

### Parallelising runs

You can run rules in parallel by providing Snakemake with multiple cores: `snakemake --cores 4 ...`.

### Debugging

You have a few options for debugging:

1. Set a `breakpoint()` inside your Python scripts and then run `snakemake ...`.
   The process will break at the `breakpoint()` and then you can use [ipdb](https://pypi.org/project/ipdb/) to inspect the current state of the variables.
1. You can set up logging inside your script.
   For logging, you can use the logger (see `scripts/dummy_script.py` for an example) and point it to a Snakemake logfile.
   Then, the rule just needs to define a logfile using the `log` option.
1. You can set up assertions in your script, which will fail when something looks wrong.
   E.g., if you always expect your data to be positive, then you can add the line `assert (df >= 0).all(), "Found negative data!"`.
   This gives you an early warning that things are wrong.
