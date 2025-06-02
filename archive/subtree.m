%subtree.m
%Generates a random spanning tree on a full subgraph of the Chimera graph
%By: Spencer Hall
%Mississippi State University and Forschungzentrum Julich

clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Preliminary Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%test to see if the library is set up correctly (check SAPI version)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define the chimera graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create a null matrix for the chimera graph (adjacency matrix)
for i = 1:n
    for j = i:n
        chimera(i,j) = 0;
    end
end

%define the chimera graph of a fully working, n-qubit chip by constructing its adjacency matrix
for row = 1:y
    for col = 1:x
        for k1 = 1:4 %vertices 1-4 of a K44 qubit group
            for k2 = 5:8 %vertices 5-8 of a K44 qubit group
                
                %define the K44 intra-connectivity
                i = k1 + (col-1)*8 + (row-1)*x*8; %vertex
                j = k2 + (col-1)*8 + (row-1)*x*8; %vertex
                chimera(i,j) = 1; %edge
                
                %define the row inter-connectivity (column intra-connectivity)
                if row ~= 8 %check for bottom edge
                    chimera(i,i+64) = 1; %edge
                end
                
                %define the column inter-connectivity (row intra-connectivity)
                if col ~= 8 %check for right edge
                    chimera(j,j+8) = 1; %edge
                end
                
            end
        end
    end
end

%inactivate the non-working qubits
for i = 1:length(qnw)
    chimera(:,qnw(i)) = 0;
    chimera(qnw(i),:) = 0;
end

%inactivate the non-working couplers
for i = 1:size(cnw,1)
    chimera(cnw(i,1),cnw(i,2)) = 0;
end

%x0 = input('X-coordinate of ulc K44 group: ')
%y0 = input('Y-coordinate of ulc K44 group: ')
ulc = input('Qubit in the upper left corner: ');
a = input('Column-dimension of rectangle: ');
b = input('Row-dimension of rectangle: ');

%ulc = 1 + 8*x*(y0-1) + 8*(x0-1);
lrc = ulc + 8*x*(b-1)+ 8*a - 1;

%inactivate qubits 1 through ulc-1
for i = 1:ulc-1
    chimera(:,i) = 0;
    chimera(i,:) = 0;
end

%inactivate qubits lrc+1 through 512
for i = lrc+1:512
    chimera(:,i) = 0;
    chimera(i,:) = 0;
end

%inactivate the rest of the qubits on the sides
for k = 1:a-1
    for i = ulc + 8*a + 64*(k-1):ulc + 8*x*k - 1
        J(i,:) = 0;
        J(:,i) = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define a random spanning tree on the chimera graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create a null matrix for the spanning tree/interactions (adjacency matrix)
for i = 1:n
    for j = i:n
       J(i,j) = 0;
    end
end

%randomly choose the starting vertex
r(1) = randi(48);
r(2) = randi([65,112]);
r(3) = randi([129,176]);
r(4) = randi([193,240]);
r(5) = randi([257,304]);
r(6) = randi([321,368]);
old = r(randi(length(r)));

%check to make sure that the starting vertex is not a non-working qubit
check = 0;
while check == 0
    for i = 1:length(qnw)
        if old == qnw(i)
            r(1) = randi(48);
            r(2) = randi([65,112]);
            r(3) = randi([129,176]);
            r(4) = randi([193,240]);
            r(5) = randi([257,304]);
            r(6) = randi([321,368]);
            old = r(randi(length(r)));
            if i == length(qnw)
                i = 0;
            end
            break
        end
        if i == length(qnw)
            check = 1;
        end
    end
end

%flag that you have visited the starting vertex
J(old,old) = 1;

%define a random spanning tree by constructing its adjacency matrix, via a random walk on the chimera graph
diagonal = diag(J);
while sum(diagonal) ~= 288-6 %check to see if you have visited all the vertices
    i = 0;
    moves = [];
    
    %look for edges to move along
    if old ~= n %check for the right boundary in the J matrix (whether or not you are at vertex n)
        for j = old+1:n %look for edges row-wise in the J matrix
            if chimera(old,j) ~= 0 %found an edge
                i = i + 1;
                moves(i) = j; %index possible moves
            end
        end
    end
    
    %look for edges to move along
    if old ~= 1 %check for the top boundary in the J matrix (whether or not you are at vertex 1)
        for j = 1:old-1 %look for edges column-wise in the J matrix
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

%check the condition for a tree ("If T is a tree, then e(T) = v(T) - 1.")
%if sum(sum(J)) ~= v-1
    %display('ERROR: This is not a tree.')
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Assign the interactions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I = input('Interaction (element of [-1.5,1], ''glass'', or ''random''): ');

%define the Edwards-Anderson glass interaction, or . . .
if strcmp(I,'glass') == 1
    gI = [-1 1];

%. . . define the bounds of the random distribution, or . . .
elseif strcmp(I,'random') == 1
    %instructions
    display('Define the bounds of the random distribution.')

    uI = input('  Upper bound (element of [-1.5,1]): ');
    %check that the upper bound is possible
    while or(uI > 1, uI < -1.5)
        display('  Upper bound must be an element of [-1.5,1].')
        uI = input('  Upper bound (element of [-1.5,1]): ');
    end

    lI = input('  Lower bound (element of [-1.5,upper bound]): ');
    %check that the lower bound is possible
    while or(lI > uI, lI < -1.5)
        display('  Lower bound must be an element of [-1.5,upper bound].')
        lI = input('  Lower bound (element of [-1.5,upper bound]): ');
    end
    
%. . . check that the given interaction is possible
else
    while or(I > 1, I < -1.5)
        display('Interaction must be an element of [-1.5,1].')
        I = input('Interaction (element of [-1.5,1]): ');
    end
end

%remove the self-interactions (the vertex flags left over from constructing the adjaceny matrix of the spanning tree)
for i = 1:n
    J(i,i) = 0;
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
save J.mat J

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