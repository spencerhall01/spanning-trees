%chimeradata.m
%Imports data from a file, formats the data, and submits the problem
%By: Spencer Hall
%Mississippi State University and Forschungzentrum Julich

clear all
clc

pause on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Preliminary Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%test to see if the library is set up correctly
display(['SAPI version: ', sapiVersion])

%test the url
%url = input('SAPI url: ');
url = 'https://qubist.dwavesys.com/sapi/';
%urlread(url);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set up a connection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create a remote SAPI connection handle
%token = input('API token: ');
token = 'your_api_token_goes_here'; %SECURITY: Do not commit real API tokens to archive - jsh 5/31/2025
%conn = sapiRemoteConnection(url, token);

%list the SAPI connection's available solvers
%sapiListSolvers(conn)

%create a SAPI solver handle
%solvername = input('Solver name: ');
solvername = 'SR10-V6';
%solver = sapiSolver(conn, solvername);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define the topology of the chip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%properties of the connected solver
%sapiSolverProperties(solver)

x = input('Number of rows of K44 arrangements: ');
y = input('Number of columns of K44 arrangements: ');

%number of qubits
n = 8*x*y;

%qnw = input('Non-working qubits, counting from 0 (vector): ');
qnw = [175,176,180,213,247,249,253,267,271,307,308,354,356,430,431,440]; %for SR10-V6 (counting from 0)
qnw = qnw + 1; %count from 1

%cnw = input('Non-working couplers (): ');
cnw(1,1) = 76; %for SR10-V6 (counting from 0)
cnw(1,2) = 84; %for SR10-V6 (counting from 0)
cnw = cnw + 1; %count from 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Specify the run-time parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Note: Default answer mode is 'histogram'.
%Note: Default number of reads is 10.
%Note: Default maximum number of answers is 1 000 for 'histogram' and num_reads for all answers ('raw').
%Note: Default annealing time is 20 us.
%Note: Default post programming thermalization time is 1 000 us.
%Note: Default post readout thermalization time is 0 us.
%Note: Default automatic scaling is true.

%answer_mode = input('Answer mode (''histogram'' or ''raw''): ');
answer_mode = 'histogram';
%num_reads = input('Number of reads (element of [1,1 000]): ');
num_reads = 1000;
if strcmp(answer_mode,'histogram') == 1
    %max_answers = input('Maximum number of answers (element of [1,1 000]): ');
    max_answers = 1000;
else
    max_answers = input('Maximum number of answers (element of [1,number of reads]): ');
end

%display('NOTE: Maximum job duration excluding readout time is 1 000 000.0 us.')
%annealing_time = input('Annealing time per read (element of [1,20 000] in us): ');
annealing_time = 20;
%programming_thermalization = input('Post programming thermalization time per read (element of [0,10 000] in us): ');
programming_thermalization = 1000;
%readout_thermalization = input('Post readout thermalization time per read (element of [0,10 000] in us): ');
readout_thermalization = 0;

%auto_scale = input('Automatic scaling (true or false): ');
auto_scale = false;

%store the parameters, including the date and time
parameters = [datenum(clock); num_reads; annealing_time; programming_thermalization; readout_thermalization; auto_scale;]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BEGIN PROGRAM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%import data from a file
uiimport -file

%wait for the user to press any key
pause %(Note: THIS COMMAND IS NECESSARY AT THIS POINT!)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Format the data for submission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%assign the magnetic field
for i = 1:n
    h(i) = 0;
end

%make sure there are no self-interaction entries
for i = 1:n
    J(i,i) = 0;
end

%make sure the adjacency matrix is upper-triangular
for i = 1:n
    for j = 1:n
        if J(i,j) ~= 0
            J(j,i) = 0;
            J(i,j) = -1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Submit the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solve the Ising problem, passing in the parameters
%answer = sapiSolveIsing(solver, h, J, 'programming_thermalization', programming_thermalization, 'num_reads', num_reads, 'max_answers', max_answers, 'readout_thermalization', readout_thermalization, 'answer_mode', answer_mode, 'auto_scale', auto_scale, 'annealing_time', annealing_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Format and save the returned answer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%extract the returned data from answer
%solutions = getfield(answer,'solutions');
%energies = getfield(answer, 'energies');
%num_occurrences = getfield(answer,'num_occurrences');
%timing = getfield(answer,'timing');

%store the lowest-energy solution
%gs = solutions(:,1);

%find and store the errors in the lowest-energy solution
j = 0;
%for i = 1:512
%    if solutions(i,1) ~= 1
%        j = j + 1;
%        errors(j,1) = i;
%        errors(j,2) = solutions(i,1);
%    end
%end

%format the returned data
%energies = energies';
%num_occurrences = num_occurrences';
%timing = timing';

%store the lowest energy, its number of occurrences, and the date & time
%lowest(1,1) = energies(1);
%lowest(1,2) = num_occurrences(1);

%check to see if the there are two solutions associated with the lowest energy
%if energies(1) == energies (2)
%    output(1,2) = num_occurrences(1) + num_occurrences(2);
%else
%    output(1,2) = num_occurrences(1);
%end

%save the solver's name, the run-time parameters, and the returned data to files
%save output.dat solvername parameters energies num_occurrences timing -ASCII
%save lowest.dat output -ASCII