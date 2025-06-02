clear all
clc

%software version
sapiVersion

%create a local SAPI connection handle
%(use to connect to local simulators)
%conn = sapiLocalConnection;

%create a remote SAPI connection handle
url = 'https://qfe.nas.nasa.gov/sapi/';
token = 'abc123';
conn = sapiRemoteConnection(url, token);

%specify the solver to use
solvername = 'C12';

%create a SAPI solver handle
solver = sapiSolver(conn, solvername);

%properties of the connected solver
props = sapiSolverProperties(solver);

%extract the properties from props
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

x = 12; %number of columns on chip
y = 12; %number of rows on chip
n = 8*x*y; %number of qubits

%used to create a rectangular subgragh
%(the example numbers given create an 8x8 graph on the chip which has its
%upper-left corner located at the first qubit)
ulc = 1; %qubit at upper-left corner of rectangle (count from 1)
a = 8; %column-dimension of rectangle
b = 8; %row-dimension of rectangle
lrc = ulc + 8*x*(b-1)+ (8*a - 1); %qubit at lower-right corner of rectangle

%non-working qubits
qnw = [26; 49; 52; 53; 98; 103; 105; 110; 168; 175; 201; 208; 215; 304; 338; 341; 372; 402; 405; 410; 413; 447; 455; 480; 485; 534; 539; 540; 576; 599; 602; 606; 642; 647; 717; 833; 853; 854; 883; 895; 905; 911; 936; 938; 943; 997; 1006; 1014; 1033; 1034; 1036; 1042; 1047; 1094; 1096]; %(counting from 0) for solver C12
qnw = qnw + 1; %count from 1

%non-working couplers
cnw =[344,351]; %(counting from 0) for solver C12
cnw = cnw + 1; %count from 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define the chimera graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create a null matrix for the graph
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
%Create J matrix (example)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J = zeros(n,n); %preallocate
for i = 1:n
    for j = i+1:n
        if chimera(i,j) == 1
            J(i,j) = -1;
        end
    end
end

%inactivate couplers associated with the non-working qubits
if isempty(qnw) ~= 1
    for i = 1:length(qnw)
        J(:,qnw(i)) = 0;
        J(qnw(i),:) = 0;
    end
end

%inactivate the non-working couplers
if isempty(cnw) ~= 1
    for i = 1:size(cnw,1)
        J(cnw(i,1),cnw(i,2)) = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create h vector (example)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = zeros(n,1); %preallocate

%inactivate the non-working qubits
if isempty(qnw) ~= 1
    for i = 1:length(qnw)
        h(qnw(i)) = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Specify the run-time parameters
%SOME INFORMATION IN THIS SECTION MAY BE OUT OF DATE!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Note: Default answer mode is 'histogram'.
%Note: Default number of reads is 10.
%Note: Default maximum number of answers is 1 000 for 'histogram' and num_reads for all answers ('raw').
%Note: Default annealing time is 20 us.
%Note: Default post programming thermalization time is 1 000 us.
%Note: Default post readout thermalization time is 0 us.
%Note: Default automatic scaling is true.

answer_mode = 'histogram'; %histogram or raw
num_reads = 1000; %element of [1,1000]
if strcmp(answer_mode,'histogram') == 1
    max_answers = 1000; %element of [1,1000]
else
    max_answers = input('Maximum number of answers (element of [1,num_reads]): ');
end

%NOTE: Maximum job duration excluding readout time is 1,000,000.0 us
annealing_time = 20; %element of [1,20000]
programming_thermalization = 1000; %element of [0,10000]
readout_thermalization = 0; %element of [0,10000] in us
auto_scale = false; %true or false

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Submit the problem and extract the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solve the Ising problem, passing in the parameters
answer = sapiSolveIsing(solver, h, J, 'programming_thermalization', programming_thermalization, 'num_reads', num_reads, 'max_answers', max_answers, 'readout_thermalization', readout_thermalization, 'answer_mode', answer_mode, 'auto_scale', auto_scale, 'annealing_time', annealing_time);

%extract the returned data from 'answer'
solutions = answer.solutions;
energies = answer.energies;
num_occurrences = answer.num_occurrences;
timing = answer.timing;

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