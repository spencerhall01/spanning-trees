%mirror.m
%Inserts a mirror into a D-Wave n-qubit chip with defined Ising parameters
%By: Spencer Hall
%Mississippi State University and Forschungzentrum Julich

clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Preliminary Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%test to see if the library and MATLAB path are set up correctly
display(['SAPI version: ', sapiVersion])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set up a connection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

conntype = input('Local connection (''l'') or remote connection (''r'')? ');
%force an acceptable input for the variable conntype
while and(strcmp(conntype,'l')~=1, strcmp(conntype,'r')~=1)
    conntype = input(' You must type ''l'' to set up a local connection or ''r'' to set up a remote connection: ');
end

%create a SAPI connection handle
if conntype == 'l'
    %create a local SAPI connection handle
    conn = sapiLocalConnection;
else
    %test url
    url = input('SAPI url: ');
    %url = 'https://qubist.dwavesys.com/sapi/';
    urlread(url);
    
    %create a remote SAPI connection handle
	token = input('API token: ');
	token = 'your_api_token_goes_here'; %SECURITY: Do not commit real API tokens to archive - jsh 5/31/2025
    conn = sapiRemoteConnection(url, token);
end

%list the SAPI connection's available solvers
solverNames = sapiListSolvers(conn);
display('Available solvers: ')
disp(solverNames)

solvername = input('Solver name: ');
%force an acceptable input for the variable solvername
flag = 0;
while flag == 0
    for i = 1:length(solverNames)
        if strcmp(solvername,solverNames(i)) == 1
            flag = 1;
        end
    end
    if flag == 0
        solvername = input(' The solver name must match one of the available solvers: ');
    end
end

%create a SAPI solver handle
solver = sapiSolver(conn, solvername);

%create a BlackBoxSolver
%blackBoxSolver = sapiBlackBoxSolver(solver);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define the topology of the chip and the section on which to compute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%properties of the connected solver
props = sapiSolverProperties(solver);
display('Solver properties: ')
disp(props)

%properties of the BlackBoxSolver
%blackBoxSolverProperties = sapiSolverProperties(blackBoxSolver);

%extract the properties
couplers = getfield(props,'couplers');
default_beta = getfield(props,'default_beta');
num_qubits = getfield(props,'num_qubits');
qubits = getfield(props,'qubits');
supported_problem_types = getfield(props,'supported_problem_types');
var_order = getfield(props,'var_order');

%store the properties
properties = [couplers, default_beta, num_qubits, qubits, supported_problem_types, var_order];

%size of the chip
x = input('Number of columns of K44 arrangements: ');
y = input('Number of rows of K44 arrangements: ');

%number of qubits
n = 8*x*y;

%non-working qubits
%qnw = input('Non-working qubits, counting from 0 (columnn vector): ');
qnw = [26; 49; 52; 53; 98; 103; 105; 110; 168; 175; 201; 208; 215; 304; 338; 341; 372; 402; 405; 410; 413; 447; 455; 480; 485; 534; 539; 540; 576; 599; 602; 606; 642; 647; 717; 833; 853; 854; 883; 895; 905; 911; 936; 938; 943; 997; 1006; 1014; 1033; 1034; 1036; 1042; 1047; 1094; 1096]; %(counting from 0) for C12
qnw = qnw + 1; %count from 1

%non-working couplers
cnw = input('Non-working couplers (two-column matrix): ');
%cnw = []; %(counting from 0)
cnw = cnw + 1; %count from 1

display('Define the rectangle on which to compute.')
%x0 = input(' Column-coordinate of upper-left-corner K44 group: ')
%y0 = input(' Row-coordinate of upper-left-corner K44 group: ')
ulc = input(' Qubit in the upper-left corner of the rectangle: ');
a = input(' Column-dimension of the rectangle: ');
b = input(' Row-dimension of the rectangle: ');

%ulc = 1 + 8*x*(y0-1) + 8*(x0-1);
lrc = ulc + 8*x*(b-1)+ (8*a - 1);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define the graph and assign the interactions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create a null matrix for the graph/interactions (adjacency matrix)
for i = 1:n
    for j = i:n
        J(i,j) = 0;
    end
end

I = input('Interaction (element of [-1.5,1], ''glass'', or ''random''): ');
    
%define the Edwards-Anderson glass interaction, or . . .
if strcmp(I,'glass') == 1
    gI = [-1 1];

%. . . define the bounds of the random distribution, or . . .
elseif strcmp(I,'random') == 1
    %instructions
    display(' Define the bounds of the random distribution.')

    lI = input('  Lower bound (element of [-1.5,1]): ');
    %check that the lower bound is possible
    while or(lI < -1.5, lI > 1)
        lI = input('   Lower bound must be an element of [-1.5,1]): ');
    end

    uI = input('  Upper bound (element of [lower bound,1]): ');
    %check that the upper bound is possible
    while or(uI < lI, uI > 1)
        uI = input('   Upper bound must be an element of [lower bound,1]): ');
    end

%. . . check that the given interaction is possible
else
    while or(I < -1.5, I > 1)
        I = input(' Interaction must be an element of [-1.5,1]): ');
    end
end

%define the chimera graph of a fully working, n-qubit chip by constructing its adjacency matrix and assign the interactions
for row = 1:y
    for col = 1:x
        for k1 = 1:4 %red vertices in a K44 arrangement
            for k2 = 5:8 %blue vertices in a K44 arrangement
                
                %define the intra-connectivity of each K44 arrangement (bipartite, fully connected) and assign the interactions
                i = k1 + 8*(col-1) + 8*x*(row-1); %red vertex/qubit
                j = k2 + 8*(col-1) + 8*x*(row-1); %blue vertex/qubit
                if strcmp(I,'glass') == 1
                    J(i,j) = gI(randi(2)); %edge/interaction
                elseif strcmp(I,'random') == 1
                    J(i,j) = lI + (uI-lI)*rand; %edge/interaction
                else
                    J(i,j) = I; %edge/interaction
                end
                
                %define the row inter-connectivity (column intra-connectivity) and assign the interactions
                if row ~= y %check for the bottom boundary
                    if strcmp(I,'glass') == 1
                        J(i,i+8*x) = gI(randi(2)); %edge/interaction
                    elseif strcmp(I,'random') == 1
                        J(i,i+8*x) = lI + (uI-lI)*rand; %edge/interaction
                    else
                        J(i,i+8*x) = I; %edge/interaction
                    end
                end
                
                %define the column inter-connectivity (row intra-connectivity) and assign the interactions
                if col ~= x %check for the right boundary
                    if strcmp(I,'glass') == 1
                        J(j,j+8) = gI(randi(2)); %edge/interaction
                    elseif strcmp(I,'random') == 1
                        J(j,j+8) = lI + (uI-lI)*rand; %edge/interaction
                    else
                        J(j,j+8) = I; %edge/interaction
                    end
                end
                
            end
        end
    end
end

%inactivate the non-working qubits
if length(qnw) ~= 0
    for i = 1:length(qnw)
        J(:,qnw(i)) = 0;
        J(qnw(i),:) = 0;
    end
end

%inactivate the non-working couplers
if length(cnw) ~= 0
    for i = 1:size(cnw,1)
        J(cnw(i,1),cnw(i,2)) = 0;
    end
end

%inactivate qubits 1 through ulc-1
if ulc ~= 1
    for i = 1:ulc-1
        J(:,i) = 0;
        J(i,:) = 0;
    end
end

%inactivate qubits lrc+1 through n
if lrc ~= n
    for i = lrc+1:n
        J(:,i) = 0;
        J(i,:) = 0;
    end
end

%inactivate the rest of the qubits on the left and right sides of the rectangle
if a ~= x
    if b ~= 1
        for k = 1:b-1
            for i = ulc + 8*a + 8*x*(k-1):ulc + 8*x*k - 1
                J(i,:) = 0;
                J(:,i) = 0;
            end
        end
    end
end

%save the interactions to a file
save J.mat J -ASCII

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Assign the magnetic fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = input('Magnetic field (element of [-2,2] or type ''random''): ');

%check that the given magnetic field is possible
while or(H > 2, H < -2)
    if H == 'random'
        break
    end
    H = input(' Magnetic field must be an element of [-2,2] or type ''random''): ');
end

%input the bounds for the random distribution
if H == 'random'
    %instructions
    display(' Define the bounds of the random distribution.')

    lH = input('  Lower bound (element of [-2,2]): ');
    %check that the lower bound is possible
    while or(lH < -2, lH > 2)
        lH = input('   Lower bound must be an element of [-2,upper bound]): ');
    end

    uH = input('  Upper bound (element of [-2,2]): ');
    %check that the upper bound is possible
    while or(uH < lH, uH > 2)
        uH = input('   Upper bound must be an element of element of [lower bound,2]): ');
    end
end

%assign the magnetic fields
if H == 'random'
    h = lH + (uH-lH).*rand(n,1);
else
    for i = 1:n
        h(i) = H;
    end
end

%inactivate non-working qubits
if length(qnw) ~= 0
    for i = 1:length(qnw)
        h(qnw(i)) = 0;
    end
end

%save the magnetic fields to a file
save h.dat h -ASCII

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Insert the vertical mirror (horizontal reflection)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create a null vector for the images
for i = 1:n
    image(i) = 0;
end

%reflect the red vertices/qubits
for col = 1:x/2
    for row = 1:y
        for k1 = 1:4
            i = k1 + 8*(col-1) + 8*x*(row-1); %red vertex/qubit
            iprime = ((5-k1) + 8*x*(row-1)) + 8*(x-1) - 8*(col-1); %red image
            image(i) = iprime; %record image in image vector
            h(iprime) = h(i); %reflect
        end
    end
end

%reflect the blue vertices/qubits
for col = 1:x/2
    for row = 1:y
        for k2 = 5:8
            j = k2 + 8*(col-1) + 8*x*(row-1); %blue vertex/qubit
            jprime = (k2 + 8*x*(row-1)) + 8*(x-1) - 8*(col-1); %blue image
            image(j) = jprime; %record image in image vector
            h(jprime) = h(j); %reflect
        end
    end
end

%fill out the rest of the image vector
for i = 1:n
    if image(i) == 0
        for j = 1:n
            if image(j) == i
                image(i) = j;
            end
        end
    end
end

%reflect the edges/interactions
for i = 1:n
    for j = i:n
        if J(i,j) ~= 0
            J(image(i),image(j)) = J(i,j); %reflect
        end
    end
end         

%save the images to a file
save Jimage.dat J -ASCII
save himage.dat h -ASCII

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END PROGRAM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Submit the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solve the Ising problem, passing in the parameters
answer = sapiSolveIsing(solver, h, J, 'programming_thermalization', programming_thermalization, 'num_reads', num_reads, 'max_answers', max_answers, 'readout_thermalization', readout_thermalization, 'answer_mode', answer_mode, 'auto_scale', auto_scale, 'annealing_time', annealing_time);

%Solve the Ising problem asynchronously
%result = sapiAsyncSolveIsing(solver, h, J);
%while ~sapiAsyncDone(result)
%    pause(1);
%end
%answer = sapiAsyncResult(result);

%solve a QUBO problem, passing in the parameters
%answer = sapiSolveQubo(solver, Q);

%Solve a QUBO problem asynchronously
%result = sapiAsyncSolveQubo(solver, Q);
%while ~sapiAsyncDone(result)
%    pause(1);
%end
%answer = sapiAsyncResult(result);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Format and save the returned answer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%extract the returned data from answer
solutions = getfield(answer,'solutions');
energies = getfield(answer,'energies');
num_occurrences = getfield(answer,'num_occurrences');
timing = getfield(answer,'timing');

%store the lowest-energy solution
gs = solutions(:,1);

%find and store the errors in the lowest-energy solution
j = 0;
for i = 1:n
    if solutions(i,1) ~= 1
        j = j + 1;
        errors(j,1) = i;
        errors(j,2) = solutions(i,1);
    end
end

%format the returned data
energies = energies';
num_occurrences = num_occurrences';
timing = timing';

%store the lowest energy, its number of occurrences, and the date & time
lowest(1,1) = energies(1);
lowest(1,2) = num_occurrences(1);
lowest(1,3) = datenum(clock);

%check to see if the there are two solutions associated with the lowest energy
if length(energies) ~= 1
    if energies(1) == energies (2)
        lowest(1,2) = num_occurrences(1) + num_occurrences(2);
    else
        lowest(1,2) = num_occurrences(1);
    end
end

%save the solver's name and properties, the run-time parameters, and the returned data to files
%save output.dat solvername parameters energies solutions num_occurrences timing -ASCII
save output.dat solvername parameters solutions energies num_occurrences -ASCII
save properties.dat properties
save lowest.dat lowest -ASCII

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot the Graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create a null matrix for the graph (adjacency matrix)
for i = 1:n
    for j = i:n
        adj(i,j) = 0;
    end
end

%define the graph by constructing its adjacency matrix
for i = 1:n
    adj(i,i) = 1; %include reflexive loops so that every vertex/qubit is plotted
    for j = i:n
        if J(i,j) ~= 0
            adj(i,j) = 1;
        end
    end
end

%inactivate the non-working qubits
for i = 1:length(qnw)
    adj(qnw(i),qnw(i)) = 0;
end

%associate each vertex/qubit with a Cartesian coordinate
for row = 1:y
    for col = 1:x
        for k1 = 1:4
            i = k1 + 8*(col-1) + 8*x*(row-1);
            j = (k1+4) + 8*(col-1) + 8*x*(row-1);
            coord(i,1) = col+2*(col-1);
            coord(i,2) = k1+5*(row-1);
            coord(j,1) = col+2*(col-1)+1;
            coord(j,2) = k1+5*(row-1);
        end
    end
end

%plot the graph
gplot(adj,coord)