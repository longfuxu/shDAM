import os
import subprocess

def run_script(script_name):
    """
    Run the specified Python script using subprocess.
    Args:
        script_name (str): Name of the script to run (e.g., 'OTdataAnalyzer.py').
    """
    try:
        print(f"Running {script_name}...")
        subprocess.run(['python', script_name], check=True)
        print(f"Successfully completed {script_name}.")
    except subprocess.CalledProcessError as e:
        print(f"Error running {script_name}: {e}")
        raise

def main():
    """
    Main function to execute OTdataAnalyzer.py, KymographAnalyzer.py, and SegementsAnalyzer.py in sequence.
    """
    # Create results directory if it doesn't exist
    results_dir = '../results'
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
        print(f"Created directory: {results_dir}")

    # Define script names in the order of execution
    scripts = ['OTdataAnalyzer.py', 'KymographAnalyzer.py', 'SegementsAnalyzer.py']
    
    # Ensure scripts exist in the current directory
    for script in scripts:
        if not os.path.isfile(script):
            raise FileNotFoundError(f"Script {script} not found in the current directory.")
    
    # Run scripts in sequence
    for script in scripts:
        run_script(script)

if __name__ == "__main__":
    print("Starting the analysis pipeline...")
    main()
    print("Analysis pipeline completed successfully. Please check the results folder.")