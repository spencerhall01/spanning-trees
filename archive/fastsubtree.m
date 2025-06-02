clear all
clc

url = 'https://qubist.dwavesys.com/sapi/';
token = 'your_api_token_goes_here'; %SECURITY: Do not commit real API tokens to archive - jsh 5/31/2025
conn = sapiRemoteConnection(url, token);
solver = sapiSolver(conn, 'SR10-V6');

x = 8;
y = 8;
n = 8*x*y;

qnw = [175,176,180,213,247,249,253,267,271,307,308,354,356,430,431,440] + 1;

cnw(1,1) = 76 + 1;
cnw(1,2) = 84 + 1;

answer_mode = 'histogram';
num_reads = 1000;
max_answers = 1000;
annealing_time = 20;
programming_thermalization = 1000;
readout_thermalization = 0;
auto_scale = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define the Chimera graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create a null matrix for the chimera graph
for i = 1:n
    for j = i:n
        chimera(i,j) = 0;
    end
end

%define the full chimera graph
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

%inactivate non-working qubits
for i = 1:length(qnw)
    chimera(:,qnw(i)) = 0;
    chimera(qnw(i),:) = 0;
end

%inactivate non-working couplers
for i = 1:size(cnw,1)
    chimera(cnw(i,1),cnw(i,2)) = 0;
end

%inactivate K44 groups 11 through 64
for i = 81:512
    chimera(:,i) = 0;
    chimera(i,:) = 0;
end

%inactivate the rest of row 1
for i = 17:64
    chimera(:,i) = 0;
    chimera(i,:) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Assign the magnetic fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:512
    h(i) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define a random spanning tree on the Chimera graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for T = 1:100

clear('J','answer','energies','num_occurrences','output')

%create a null matrix for the spanning tree
for i = 1:n
    for j = i:n
       J(i,j) = 0;
    end
end

%randomly choose the starting vertex
r(1) = randi(16);
r(2) = randi([65,80]);
old = r(randi(length(r)));

%check to make sure that the starting vertex is not a non-working qubit
check = 0;
while check == 0
    for i = 1:length(qnw)
        if old == qnw(i)
            r(1) = randi(16);
            r(2) = randi([65,80]);
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

%define a random spanning tree
diagonal = diag(J);
while sum(diagonal) ~= 32 %check to see if you have visited all the vertices
    i = 0;
    moves = [];
    
    if old ~= n %check for the right boundary in the J matrix (whether or not you are at vertex n)
        for j = old+1:n %look for edges row-wise in the J matrix
            if chimera(old,j) ~= 0 %found an edge
                i = i + 1;
                moves(i) = j; %index possible moves
            end
        end
    end
    
    if old ~= 1 %check for the top boundary in the J matrix (whether or not you are at vertex 1)
        for j = 1:old-1 %look for edges column-wise in the J matrix
            if chimera(j,old) ~= 0 %found an edge
                i = i + 1;
                moves(i) = j; %index posible moves
            end
        end
    end
    
    new = moves(randi(length(moves))); %next vertex
    
    if diagonal(new) == 0 %check to see if you have already visited the next vertex
        J(min(old,new),max(old,new)) = 1; %flag the edge to move along
        J(new,new) = 1; %flag the vertices that you have visited
        diagonal = diag(J);
    end
    
    old = new; %move
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Assign the interactions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%remove self-interactions
for i = 1:n
    J(i,i) = 0;
end

%Assign the interactions
gI = [-1 1];
for i = 1:n
    for j = i:n
        if J(i,j) == 1
            J(i,j) = gI(randi(2));
        end
    end
end

%save the tree
save(['EAtrees2x2_',int2str(T)], 'J')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Submit the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear('i','j','old','check','diagonal','moves','new','interaction','R')

for R = 1:100

clear('answer','energies','num_occurrences','output')
clc
[T R]

%solve the Ising problem, passing in the parameters
answer = sapiSolveIsing(solver, h, J, 'programming_thermalization', programming_thermalization, 'num_reads', num_reads, 'max_answers', max_answers, 'readout_thermalization', readout_thermalization, 'answer_mode', answer_mode, 'auto_scale', auto_scale, 'annealing_time', annealing_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Format and save the returned answer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%extract the returned data from answer
%solutions = getfield(answer,'solutions');
energies = getfield(answer, 'energies');
num_occurrences = getfield(answer,'num_occurrences');
%timing = getfield(answer,'timing');

%store the lowest-energy solution
%gs = solutions(:,1);

%find and store the errors in the lowest-energy solution
%j = 0;
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
output(1,1) = energies(1);
%output(T,2) = num_occurrences(1);
%output(3) = datenum(clock);

if energies(1) == energies (2)
    output(1,2) = num_occurrences(1) + num_occurrences(2);
else
    output(1,2) = num_occurrences(1);
end

%save the returned data to files
%save energies.dat energies -ASCII
%save num_occurrences.dat num_occurrences -ASCII
save EA100x100randspan2x2.dat output -append -ASCII

end
end