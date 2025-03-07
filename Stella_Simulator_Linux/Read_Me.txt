To Install:

1. Download the zip archive for the operating system you are using. Unzip this into a directory you want to run stella simulator from (it can be included in your path, or you can use the path when calling Stella Simulator). 

2. Download the license file (license.xml) from the my product page, and put it in the same directory as the program.

To Run: 

The command line options match those of the main program (https://iseesystems.com/help/default.htm#cshid=1195) except -r is not needed, -rs is the same as -r, and -nq, -rd, and -s are all ignored.  There are also additional options which will appear if you just type the program name and press enter:

Usage: stella_simulator <options> model_file
  -0 arg               Set specified variable to 0 before first run
  -1 arg               Set specified variable to 1 before first run
  -d arg               Seconds to delay between successive runs (default: 0)
  -h arg               Name of file to use for handshakes (default: none)
  -i                   Import now before each run
  -ld arg              Load data from the file as a run
  -ltm arg             Run LTM and put results in file (and file_defs.txt)
  -ltmn arg            Exhaustive loop search threshold
  -p arg               Pause at the given interval (0: DT; default: no pause)
  -pd arg              Seconds to delay before resuming after a pause
                        (default: 0)
  -pe arg              Off|Model|Interface Observe pause events from the
                        mdodel or interface (default off)
  -ph arg              Name of file to use for pause/resume handshakes
                        (changes -pd to maximum wait/overrides -pf; default:
                        none)
  -pi                  Pause immediately after initialization
  -pir                 Pause immediately after initialization, then
                        reinitialize on resume
  -ps arg              Name of the status file (default: none). Use with -ph
  -q                   Quiet mode (only errors are output)
  -r                   Run the model once (default)
  -rn arg              Run specified number of times
  -ro                  Perform optimization (requires optimization be set up
                        in the model)
  -rs                  Perform sensitivity (requires sensitivity be set up in
                        the model)
  -x                   Export now after each run

These options are there to allow you to control how often it exports and when it resumes after a pause (waiting for you to update an import sheet).  For example:

   stella_simulator -p 4 -pd 0.25 my_model.stmx

will pause every 4 time units, export, delay 1/4 second, import, and then resume.

Handshakes:
- Stella Simulator erases the handshake file after the export is complete
- The client must write some data (any data - I just write the character '1') to this same file after it has updated the import sheet
- You can abort the simulation by writing **#** to the handshake file
