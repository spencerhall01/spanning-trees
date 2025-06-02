%WARNING: THIS CODE DOES NOT WORK AS INTENDED . . . ARCHIVED FOR LATER REFERENCE. - jsh 5/31/2025

%quadraticforms.m
%Defines a graph on a D-Wave n-qubit chip and assigns Ising parameters
%By: Spencer Hall
%Mississippi State University and Forschungzentrum Julich

clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Preliminary Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%test to see if the library and MATLAB path are set up correctly
display(['SAPI version: ', sapiVersion])

%test url
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
%Define the topology of the chip and the section on which to compute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%properties of the connected solver
%sapiSolverProperties(solver)

x = input('Number of columns of K44 arrangements: ');
y = input('Number of rows of K44 arrangements: ');

%number of qubits
n = 8*x*y;

%qnw = input('Non-working qubits, counting from 0 (columnn vector): ');
qnw = []; %for simulator (none) (counting from 0)
qnw = qnw + 1; %count from 1

%cnw = input('Non-working couplers (two-column matrix): ');
cnw(1,1) = 76; %for SR10-V6 (counting from 0)
cnw(1,2) = 84; %for SR10-V6 (counting from 0)
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
%Assign the ground state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GS = [-1 1];
for i = 1:512
    gs(i) = GS(randi(2));
end

%inactivate non-working qubits
for i = 1:length(qnw)
    gs(qnw(i)) = 0;
end

%save the ground state to a file
save gs.dat gs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Assign the magnetic fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = input('Magnetic field (element of [-2,2] or type ''random''): ');

%check that the given magnetic field is possible
while or(H > 2, H < -2)
    if H == 'random'
        break
    end
    display('Magnetic field must be an element of [-2,2].')
    H = input('Magnetic field (element of [-2,2] or type ''random''): ');
end

%input the bounds for the random distribution
if H == 'random'
    %instructions
    display('Define the bounds of the random distribution.')
    
    uH = input('  Upper bound (element of [-2,2]): ');
    %check that the upper bound is possible
    while or(uH > 2, uH < -2)
        display('  Upper bound must be an element of [-2,2].')
        uH = input('  Upper bound (element of [-2,2]): ');
    end
    
    lH = input('  Lower bound (element of [-2,upper bound]): ');
    %check that the lower bound is possible
    while or(lH > uH, lH < -2)
        display('  Lower bound must be an element of [-2,upper bound].')
        lH = input('  Lower bound (element of [-2,upper bound]): ');
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
for i = 1:length(qnw)
    h(qnw(i)) = 0;
end

%save the magnetic fields to a file
save h.dat h

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define the graph and assign the interactions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create a null matrix for the graph/interactions (adjacency matrix)
for i = 1:n
    for j = i:n
        J(i,j) = 0;
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
                if gs(i)==gs(j)
                    J(i,j) = -1; %edge/interaction
                else
                    J(i,j) = 1; %edge/interaction
                end
                
                %define the row inter-connectivity (column intra-connectivity) and assign the interactions
                if row ~= y %check for the bottom boundary
                    if gs(i)==gs(j)
                        J(i,j) = -1; %edge/interaction
                    else
                        J(i,j) = 1; %edge/interaction
                    end
                end
                
                %define the column inter-connectivity (row intra-connectivity) and assign the interactions
                if col ~= x %check for the right boundary
                    if gs(i)==gs(j)
                        J(i,j) = -1; %edge/interaction
                    else
                        J(i,j) = 1; %edge/interaction
                    end
                end
                
            end
        end
    end
end

%inactivate the non-working qubits
for i = 1:length(qnw)
    J(:,qnw(i)) = 0;
    J(qnw(i),:) = 0;
end

%inactivate the non-working couplers
for i = 1:size(cnw,1)
    J(cnw(i,1),cnw(i,2)) = 0;
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
if a ~= 8
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
save J.mat J

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
%    lowest(1,2) = num_occurrences(1) + num_occurrences(2);
%else
%    lowest(1,2) = num_occurrences(1);
%end

%save the solver's name, the run-time parameters, and the returned data to files
%save output.dat solvername parameters energies num_occurrences timing -ASCII
%save lowest.dat output -ASCII

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