%randspan.m
%Defines a random spanning tree on the Chimera graph and assigns the parameters
%By: Spencer Hall
%Mississippi State University and Forschungzentrum Julich

clear all
clc

display('randspan.m')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%List of user-defined variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%conntype (Local or remote connection? l/r) %checked
%url (SAPI url)
%token (API token)
%solvername (name of solver to use) %checked
%qsage (Create QSage solver? y/n) %checked

%x (number of columns on chip)
%y (number of rows on chip)
%qnw (zero-based indices of qubits not working (vector))
%cnw (zero-based indices of couplers not working (two-column vector))
%x0 (column-coordinate of upper left corner of subgraph)
%y0 (row-coordinate of upper left corner of subgraph)
%ulc (qubit in upper left corner of subgraph)
%a (column-dimension of the subgraph)
%b (row-dimension of the subgraph)
%NOTE: {x0,y0} and {ulc} are mutually exclusive

%answer_mode
%max_answers
%num_reads
%annealing_time
%programming_thermalization
%readout_thermalization
%auto_scale

%I (interaction weights for couplers) %checked
%lI (lower bound for random weights) %checked
%uI (upper bound for random weights) %checked

%H (external magnetic field biases for qubits) %checked
%lH (lower bound for random biases) %checked
%uH (upper bound for random biases) %checked

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Test the library and MATLAB path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%test to see if the library and MATLAB path are set up correctly
display(['SAPI version: ', sapiVersion])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Connect to a solver and display its properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

conntype = input('Local or remote connection (l/r)? ', 's');
%force an acceptable input for the variable conntype
while and(strcmp(conntype,'l')~=1, strcmp(conntype,'r')~=1)
    conntype = input(' You must type l or r: ');
end

%create a SAPI connection handle
if conntype == 'l'
    %create a local SAPI connection handle
    conn = sapiLocalConnection;
else
    %input and test the url
    url = input('SAPI url: ', 's');
    %url = 'https://qfe.nas.nasa.gov/sapi/';
    %urlread(url);
    
    %create a remote SAPI connection handle
	%token = input('API token: ', 's');
	token = 'your_api_token_goes_here'; %SECURITY: Do not commit real API tokens to archive - jsh 5/31/2025
    conn = sapiRemoteConnection(url, token);
end

%list the SAPI connection's available solvers
solverNames = sapiListSolvers(conn);
display('Available solvers: ')
disp(solverNames)

solvername = input('Solver name: ', 's');
%force an acceptable input for the variable solvername
flag = 0;
while flag == 0
    for i = 1:length(solverNames)
        if strcmp(solvername,solverNames(i)) == 1
            flag = 1;
        end
    end
    if flag == 0
        solvername = input(' The solver name must match one of the available solvers: ', 's');
    end
end

%create a SAPI solver handle
solver = sapiSolver(conn, solvername);

%create a QSage solver
qsage = input('Would you like to create a QSage solver (y/n)? ', 's');
%force an acceptable input for the variable qsage
flag = 0;
while flag == 0
    if qsage == 'y'
        QSageSolver = sapiSolveQSage(solver);
        flag = 1;
    elseif qsage == 'n'
        flag = 1;
    else
        input(' You must type y or n: ');
    end
end

%properties of the connected solver
props = sapiSolverProperties(solver);
display('Solver properties: ')
disp(props)

%properties of the QSage solver
if qsage == 'y'
    qsageprops = sapiSolverProperties(QSageSolver);
    display('QSage properties: ')
    disp(qsageprops)
end

%extract the properties
annealing_time_range = props.annealing_time_range;
chip_id = props.chip_id;
couplers = props.couplers;
default_annealing_time = props.default_annealing_time;
default_programming_thermalization = props.default_programming_thermalization;
default_readout_thermalization = props.default_readout_thermalization;
h_range = props.h_range;
j_range = props.j_range;
num_qubits = props.num_qubits;
parameters = props.parameters;
programming_thermalization_range = props.programming_thermalization_range;
qubits = props.qubits;
server_version = props.server_version;
supported_problem_types = props.supported_problem_types;

%store the properties
%properties = [annealing_time_range, chip_id, couplers, default_annealing_time, default_programming_thermalization, default_readout_thermalization, h_range, j_range, num_qubits, parameters, programming_thermalization_range, qubits, server_version, supported_problem_types];

%save the solver's properties to a file
save properties.dat props

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define the topology of the chip and the section on which to compute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%size of the chip
x = input('Number of columns of K44 arrangements: ');
y = input('Number of rows of K44 arrangements: ');

%number of qubits
n = 8*x*y;

%non-working qubits
%qnw = input('Non-working qubits, counting from 0 (vector): ');
qnw = [26; 49; 52; 53; 98; 103; 105; 110; 168; 175; 201; 208; 215; 304; 338; 341; 372; 402; 405; 410; 413; 447; 455; 480; 485; 534; 539; 540; 576; 599; 602; 606; 642; 647; 717; 833; 853; 854; 883; 895; 905; 911; 936; 938; 943; 997; 1006; 1014; 1033; 1034; 1036; 1042; 1047; 1094; 1096]; %(counting from 0) for C12
qnw = qnw + 1; %count from 1

%non-working couplers
%cnw = input('Non-working couplers (two-column matrix): ');
cnw =[344,351]; %(counting from 0) for C12
cnw = cnw + 1; %count from 1

display('Define the rectangle on which to compute.')
x0 = input(' Column-coordinate of upper-left-corner K44 group: ');
y0 = input(' Row-coordinate of upper-left-corner K44 group: ');
%ulc = input(' Qubit in the upper-left corner of the rectangle: ');
a = input(' Column-dimension of the rectangle: ');
b = input(' Row-dimension of the rectangle: ');

ulc = 1 + 8*x*(y0-1) + 8*(x0-1);
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

answer_mode = input('Answer mode (histogram or raw): ', 's');
%answer_mode = 'histogram';
num_reads = input('Number of reads (element of [1,1000]): ');
%num_reads = 1000;
if strcmp(answer_mode,'histogram') == 1
    max_answers = input('Maximum number of answers (element of [1,1000]): ');
    %max_answers = 1000;
else
    max_answers = input('Maximum number of answers (element of [1,number of reads]): ');
end

display('NOTE: Maximum job duration excluding readout time is 1000000.0 us.')
annealing_time = input(' Annealing time per read (element of [1,20000] in us): ');
%annealing_time = 20;
programming_thermalization = input(' Post programming thermalization time per read (element of [0,10000] in us): ');
%programming_thermalization = 1000;
readout_thermalization = input(' Post readout thermalization time per read (element of [0,10000] in us): ');
%readout_thermalization = 0;

auto_scale = input('Automatic scaling (true or false): ');
%auto_scale = false;

%store the parameters, including the date and time
%parameters = [datenum(clock); answer_mode; num_reads; max_answers; annealing_time; programming_thermalization; readout_thermalization; auto_scale]';

%save the annealing parameters to a file
%save parameters.dat parameters -ASCII

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BEGIN PROGRAM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define the chimera graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create a null matrix for the graph (adjacency matrix)
chimera = zeros(n,n); %preallocate

%define the chimera graph of a fully working, n-qubit chip by constructing its adjacency matrix
for row = 1:y
    for col = 1:x
        for k1 = 1:4 %red vertices in a K44 arrangement
            for k2 = 5:8 %blue vertices in a K44 arrangement
                
                %define the intra-connectivity of each K44 arrangement (bipartite, fully connected)
                i = k1 + 8*(col-1) + 8*x*(row-1); %red vertex
                j = k2 + 8*(col-1) + 8*x*(row-1); %blue vertex
                chimera(i,j) = 1; %edge
                
                %define the row inter-connectivity (column intra-connectivity)
                if row ~= y %check for bottom edge
                    chimera(i,i+8*x) = 1; %edge
                end
                
                %define the column inter-connectivity (row intra-connectivity)
                if col ~= x %check for right edge
                    chimera(j,j+8) = 1; %edge
                end
                
            end
        end
    end
end

%inactivate the non-working qubits
if isempty(qnw) ~= 1
    for i = 1:length(qnw)
        chimera(:,qnw(i)) = 0;
        chimera(qnw(i),:) = 0;
    end
end

%inactivate the non-working couplers
if isempty(cnw) ~= 1
    for i = 1:size(cnw,1)
        chimera(cnw(i,1),cnw(i,2)) = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Create subgraph%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%inactivate qubits 1 through ulc-1
if ulc ~= 1
    for i = 1:ulc-1
        chimera(:,i) = 0;
        chimera(i,:) = 0;
    end
end

%inactivate qubits lrc+1 through n
if lrc ~= n
    for i = lrc+1:n
        chimera(:,i) = 0;
        chimera(i,:) = 0;
    end
end

%inactivate the rest of the qubits on the left and right sides of the rectangle
if a ~= x
    if b ~= 1
        for k = 1:b-1
            for i = ulc + 8*a + 8*x*(k-1):ulc + 8*x*k-1
                chimera(i,:) = 0;
                chimera(:,i) = 0;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define a random spanning tree on the Chimera graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for T = 1:100

%clear('J','tnw','i','j','r','old','flag','diagonal','moves','new','I','gI','lI','uI','H','gH','lH','uH','h','R')

%create a null matrix for the spanning tree/interactions (adjacency matrix)
J = zeros(n,n); %preallocate

%find the non-working qubits between the ulc and the lrc and add them to the total not working on the rectangle
tnw = 0;
for i = 1:length(qnw)
    for j = ulc:lrc
        if qnw(i) == j
            tnw = tnw + 1;
        end
    end
end

%find the non-working qubits on the left and right sides of the rectangle and remove them from the total not working on the rectangle
if a ~= x
    if b ~= 1
        for i = 1:length(qnw)
            for k = 1:b-1
                for j = ulc + 8*a + 8*x*(k-1):ulc + 8*x*k - 1
                    if qnw(i) == j
                        tnw = tnw-1;
                    end
                end
            end
        end
    end
end

%randomly choose the starting vertex
r = zeros(b,1); %preallocate
for i = 1:b
    r(i) = randi([ulc + 8*x*(i-1),ulc + 8*a-1 + 8*x*(i-1)]);
end
old = r(randi([1,length(r)]));

%check to make sure that the starting vertex is not a non-working qubit
flag = 0;
while flag == 0
    for i = 1:length(qnw)
        if old == qnw(i)
            for j = 1:b
                r(j) = randi([ulc + 8*x*(j-1),ulc + 8*a-1 + 8*x*(j-1)]);
            end
            old = r(randi([1,length(r)]));
            if i == length(qnw)
                i = 0;
            end
        end
        if i == length(qnw)
            flag = 1;
        end
    end
end

%flag that you have visited the starting vertex
J(old,old) = 1;

%define a random spanning tree by constructing its adjacency matrix, via a random walk on the chimera graph
diagonal = diag(J);
while sum(diagonal) ~= 8*a*b-tnw %check to see if you have visited all the vertices
    i = 0;
    moves = [];
    
    %look for edges to move along
    if old ~= n %check for the right boundary in the chimera matrix (whether or not you are at vertex n)
        for j = old+1:n %look for edges row-wise (along a row) in the chimera matrix (row index constant/column index changes)
            if chimera(old,j) ~= 0 %found an edge
                i = i + 1;
                moves(i) = j; %index possible moves
            end
        end
    end
    
    %still looking for edges to move along
    if old ~= 1 %check for the top boundary in the chimera matrix (whether or not you are at vertex 1)
        for j = 1:old-1 %look for edges column-wise (along a column) in the chimera matrix (row index changes/column index constant)
            if chimera(j,old) ~= 0 %found an edge
                i = i + 1;
                moves(i) = j; %index posible moves
            end
        end
    end
    
    %randomly choose the next vertex
    new = moves(randi(length(moves))); %next vertex
    
    %construct the adjacency matrix for the spanning tree
    if diagonal(new) == 0 %check to see if you have already visited the next vertex
        J(min(old,new),max(old,new)) = 1; %flag the edge to move along
        J(new,new) = 1; %flag the vertices that you have visited
        diagonal = diag(J);
    end
    
    %move
    old = new;
end

%remove the self-interactions (the vertex flags left over from constructing the adjaceny matrix of the spanning tree)
for i = 1:n
    J(i,i) = 0;
end

%check a necessary, but not sufficient, condition for a graph to be a tree ("If T is a tree, then e(T) = v(T) - 1.")
if sum(sum(J)) ~= (8*a*b-tnw)-1
    display('WARNING: This is not a tree.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Assign the interactions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I = input('Interaction (element of [-1,1], glass, or random): ', 's');
    
%define the Edwards-Anderson glass interaction, or . . .
if strcmp(I,'glass') == 1
    gI = [-1 1];

%. . . define the bounds of the random distribution, or . . .
elseif strcmp(I,'random') == 1
    %instructions
    display(' Define the bounds of the random distribution.')

    lI = input('  Lower bound (element of [-1,1]): ');
    %check that the lower bound is possible
    while or(lI < j_range(1), lI > j_range(2))
        lI = input('   Lower bound must be an element of [-1,1]): ');
    end

    uI = input('  Upper bound (element of [lower bound,1]): ');
    %check that the upper bound is possible
    while or(uI < lI, uI > j_range(2))
        uI = input('   Upper bound must be an element of [lower bound,1]): ');
    end

%. . . check that the given interaction is possible
else
    while or(I < j_range(1), I > j_range(2))
        I = input(' Interaction must be an element of [-1,1]: ');
    end
end

%assign the interactions
for i = 1:n
    for j = i:n
        if J(i,j) == 1
            if strcmp(I,'glass') == 1
                J(i,j) = gI(randi(2));
            elseif strcmp(I,'random') == 1
                J(i,j) = lI + (uI-lI)*rand;
            else
                J(i,j) = I;
            end
        end
    end
end

%save the spanning tree/interactions to a file
save J.dat J -ASCII
%save(['EAtrees8x8_',int2str(T)], 'J')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Assign the magnetic fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = input('Magnetic field (element of [-2,2], glass, random): ', 's');

%define the Edwards-Anderson glass interaction, or . . .
if strcmp(H,'glass') == 1
    gH = [-1 1];

%. . . define the bounds of the random distribution, or . . .
elseif strcmp(H,'random') == 1
    %instructions
    display(' Define the bounds of the random distribution.')

    lH = input('  Lower bound (element of [-2,2]): ');
    %check that the lower bound is possible
    while or(lI < h_range(1), lI > h_range(2))
        lH = input('   Lower bound must be an element of [-2,2]): ');
    end

    uH = input('  Upper bound (element of [lower bound,2]): ');
    %check that the upper bound is possible
    while or(uI < lI, uI > h_range(2))
        uH = input('   Upper bound must be an element of [lower bound,2]): ');
    end

%. . . check that the given field is possible
else
    while or(H < h_range(1), H > h_range(2))
        H = input(' Field must be an element of [-2,2]: ');
    end
end

%assign the magnetic fields
h = zeros(n,1); %preallocate
if strcmp(H,'glass') == 1
    for i = 1:n
        h(i) = gH(randi(2));
    end
elseif strcmp(H,'random')
    h = lH + (uH-lH).*rand(n,1);
else
    for i = 1:n
        h(i) = H;
    end
end

%inactivate non-working qubits
if isempty(qnw) ~= 1
    for i = 1:length(qnw)
        h(qnw(i)) = 0;
    end
end

%save the magnetic fields to a file
save h.dat h -ASCII

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END PROGRAM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Submit the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for R = 1:100
    
%clear('time','answer','solutions','energies','num_occurrences','timing','errors','lowest')
%clc
%[T R]

%solve the Ising problem, passing in the parameters
answer = sapiSolveIsing(solver, h, J, 'programming_thermalization', programming_thermalization, 'num_reads', num_reads, 'max_answers', max_answers, 'readout_thermalization', readout_thermalization, 'answer_mode', answer_mode, 'auto_scale', auto_scale, 'annealing_time', annealing_time);

%Solve the Ising problem asynchronously
%result = sapiAsyncSolveIsing(solver, h, J, 'programming_thermalization', programming_thermalization, 'num_reads', num_reads, 'max_answers', max_answers, 'readout_thermalization', readout_thermalization, 'answer_mode', answer_mode, 'auto_scale', auto_scale, 'annealing_time', annealing_time);
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

%record the approximate time of submission (this time recorded will be later than the actual time)
timesub = datenum(clock);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Format and save the returned answer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%extract the returned data from answer
solutions = answer.solutions;
energies = answer.energies;
num_occurrences = answer.num_occurrences;
timing = answer.timing;

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

energies_occurrences = [energies num_occurrences];

%store the lowest energy, its number of occurrences, and the date & time
lowest(1,1) = energies(1);
lowest(1,2) = num_occurrences(1);
lowest(1,3) = datenum(clock);

datetime = datenum(clock);

%check to see if the there are two solutions associated with the lowest energy
if length(energies) ~= 1
    if energies(1) == energies (2)
        lowest(1,2) = num_occurrences(1) + num_occurrences(2);
    else
        lowest(1,2) = num_occurrences(1);
    end
end

%save the solver's name, the run-time parameters, and the returned data to files
%save output.dat solvername properties parameters timesub solutions energies num_occurrences timing -ASCII
save output.dat solvername timesub solutions energies num_occurrences -ASCII
save lowest.dat lowest -ASCII
%save solutions.dat solutions -ASCII
%save EA10x100randspan8x8.dat lowest -append -ASCII

%end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot the Graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create a null matrix for the graph (adjacency matrix)
adj = zeros(n,n); %preallocate

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
coord = zeros(n,2);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save the problem in a format for cutting and pasting into the web interface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this must be the first line
first_line = [num_qubits length(qubits)+length(couplers)];
save web.dat first_line -ASCII

%format the bias data
k = 1;
flag = 0;
for i=1:n
    for j = 1:length(qnw)
        if qnw(j) == i
            flag = 1; %do not write the qnw---the website doesn't like that
        end
    end
    if flag == 0
        web_qubits(k,1)=i-1;
        web_qubits(k,2)=i-1;
        web_qubits(k,3)=h(i);
        k = k + 1;
    end
    flag = 0; %set the flag back to 0 for the next round
end
save web.dat web_qubits -append -ASCII

%format the coupler data
k=1;
for i=1:n
    for j=i:1152
        if J(i,j)~=0
            web_couplers(k,1) = i-1;
            web_couplers(k,2) = j-1;
            web_couplers(k,3) = J(i,j);
            k = k + 1;
        end
    end
end
save web.dat web_couplers -append -ASCII