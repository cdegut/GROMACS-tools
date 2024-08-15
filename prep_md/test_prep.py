import subprocess

def run_shell_command(command):
    # Start the process
    process = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    output = ""
    
    try:
        while process.poll() is not None:
            # Read output line by line
            output_line_out = process.stdout.readline()
            output_line_err = process.stderr.readline()

            # Print the output to the console
            if output_line_out:
                print(output_line_out, end="")
                output += output_line_out
            
                        # Print the output to the console
            if output_line_err:
                print(output_line_err, end="")
                output += output_line_err

            # Allow the user to input responses
         # Check if the process is still running
        user_input = input()  # Get user input
        process.stdin.write(user_input + "\n")  # Send input to the process
        process.stdin.flush()  # Ensure it gets sent immediately

        # Ensure we capture the final output
        remaining_output = process.communicate()[0]
        output += remaining_output
    
    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        process.stdout.close()
        process.stdin.close()
    
    return output

command = "gmx pdb2gmx -ignh -f U53_AF_rank1.pdb -o starting_structure.gro"
captured_output = run_shell_command(command)
print("\nCaptured Output:\n", captured_output)