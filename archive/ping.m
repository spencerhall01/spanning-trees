%ping.m
%Test the library, the url, and the connection
%By: Spencer Hall
%Mississippi State University and Forschungzentrum Julich

clear all
clc

display('ping.m')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Test the library and MATLAB path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%test to see if the library is set up correctly
display(['SAPI version: ', sapiVersion])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Connect to a solver and display its properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

conntype = input('Local or remote connection (''l''/''r'')? ');
%force an acceptable input for the variable conntype
while and(strcmp(conntype,'l')~=1, strcmp(conntype,'r')~=1)
    conntype = input(' You must type ''l'' or ''r'': ');
end

%create a SAPI connection handle
if conntype == 'l'
    %create a local SAPI connection handle
    conn = sapiLocalConnection;
else
    %input and test the url
    url = input('SAPI url: ');
    %url = 'https://qfe.nas.nasa.gov/sapi/';
    %urlread(url);
    
    %create a remote SAPI connection handle
	%token = input('API token: ');
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

%create a QSage solver
qsage = input('Would you like to create a QSage solver (''y''/''n'')? ');
if qsage == 'y'
    QSageSolver = sapiSolveQSage(solver);
elseif qsage == 'n'
    break
else
    input(' You must type ''y'' or ''n'': ');
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