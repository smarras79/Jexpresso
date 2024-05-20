# Add a new test to CI:
Follow these simple steps to add your new PROBLEM_NAME/YOUR_TEST_DIR_NAME to the CI suite:

       1. Set the following keys in problems/equations/PROBLEM_NAME/YOUR_TEST_DIR_NAME/user_inputs.jl
       	  - :output_format     => "hdf5",
	  - :output_dir        => "./test/CI-ref",
	  - :loverwrite_output => true,
	                      
       2. Run your test as usual
       
       3. cp -rf problems/equations/PROBLEM_NAME/YOUR_TEST_DIR_NAME test/CI-runs/PROBLEM_NAME/

       4. Edit user_inputs.jl to replace 
       	     :output_dir => "./test/CI-ref",
	  with
    	     :output_dir => "./test/CI-runs",
       
       5. Open test/runtests.jl and add the following line if necessary:
           `@time @testset "PROBLEM_NAME" begin include("CI-runs/PROBLEM_NAME/runtests.jl") end`
           Replace PROBLEM_NAME with the one that contains YOUR_TEST_DIR_NAME. 
           Notice that you do not need to add this new line if YOUR_TEST_DIR_NAME is already there.
        
       6. cp test/CI-runs/CompEuler/thetaTracers/Tests.jl test/CI-runs/PROBLEM_NAME/YOUR_TEST_DIR_NAME/
        
       7. edit test/CI-runs/CompEuler/thetaTracers/Tests.jl and replace 
          @testset "JEXPRESSO Examples" begin run_example("CompEuler", "thetaTracers") end
          with
          @testset "JEXPRESSO Examples" begin run_example("PROBLEM_NAME", "YOUR_TEST_DIR_NAME") end
	  
	Done. At this point the CI will run when you push the code.

# General notes:
    The `runtest.jl` file is part of a Julia package's continuous integration (CI) testing framework. 
    It uses Julia's built-in `Test` module to define and execute test sets for "Jexpresso". 
    Here’s a breakdown of how it is structured and how it contributes to the CI process:

    1. **Importing Dependencies**: At the beginning of the file, there are `using` statements for `Test` and `Jexpresso`. 
        This means the test file relies on these three packages. 
        `Test` is used for writing the test cases and `Jexpresso` is the main package being tested.

    2. **`run_example`**: 
        This function takes two arguments: `parsed_equations` and `parsed_equations_case_name`, which are the names for specific test cases or example scenarios within the Jexpresso framework. 
        The function sets up an environment, navigates to a specific directory where the test case resides, clears any existing command-line arguments, and then simulates the running of a Jexpresso application by pushing the relevant arguments and including the main Jexpresso script. 
        It wraps the execution within a `@testset`, allowing for grouped test reporting and isolation. If the script runs without throwing an error, a test within this set passes; otherwise, it captures the error, prints a portion of it, and fails the test.

    3. **Setting Up Environment Variables**: 
        Inside the `run_example` function, it adjusts the `ENV["JEXPRESSO_HOME"]` to ensure it points to the correct base directory for the Jexpresso framework.

    4. **Error Handling**:
        It includes a `try-catch` block. If the Jexpresso example runs successfully, it proceeds without interruption. 
        If there's an error, it prints out a portion of the error message and fails the current test by asserting `@test false`.

    5. **Running Test Sets for Examples**: 
        After defining `run_example`, the script defines a larger `@testset` titled "JEXPRESSO Examples". 
        Within this set, it iterates over a predefined list of example scenarios, grouped by problem names and case names. 
        Each of these scenarios is passed to the `run_example` function to be run as part of the overall test suite.

    6. **Integration with CI**: 
        This file is integrated into the CI pipeline via GitHub Actions. The CI pipeline will report success if all tests pass, or failure if any test fails, preventing merging broken code into main branches or tagging releases with failing tests.

    This structure allows for automated testing of various scenarios within the Jexpresso framework, ensuring that any new changes don't break existing functionality and that all the defined examples work as expected.

