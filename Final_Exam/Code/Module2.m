% AOBD FInal Exam Paper - Career Recommender System from LinkedIn's user's profiles 
% Module 1 -> Suggesting skills that are required for particular career position
% Created and Written by Nand Parikh
% Roll No: 1401023

clc; clearvars; close all;

%% INTIALIZAION PHASE
S = randi([1 100], 3, 500); % "Skill Space" Matrix - random
Scount = zeros(1,500);      % Matrix to count how many user has particular skill
J = zeros(3,500);           % "Job Space" Matrix
Jcount = zeros(1,500);      % Matrix to count how many user has particular job
MarkJ = zeros(1,500);       % Boolean array for job
A = randi([1 10], 3, 3);    % Transformation Matrix of 3x3 -> Note: Will always be invertible

%% DATA READER MODULE
[ID,T1] = xlsread('Data.xlsx');     % Read the cleaned Data
JobTitle = T1(2:length(T1),3);      
Skills = T1(2:length(T1),2);
MapCnt=1;JobMapCnt=1;               % Counters for no of unique skills and jobs
Names = cell(1);
map = containers.Map();             % HashMap for skillset 
JobMap = containers.Map();          % HashMap for Job

%% Main Algorithm - 1
% Clustering of similar skills and mapping to respective jobs
for k=1:length(ID)  % For each profile
    Ui = zeros(1,1);% Create empty skill set of particular user
    % Iteration through Skiils of that user
    str = char(Skills(k));
    str1 = char('');
    for i=1:length(str)
        if str(i)==','
            if map.size == 0
                map(str1) = MapCnt;
                Names(MapCnt) = cellstr(str1);
                MapCnt=MapCnt+1;
            else
                tf = isKey(map,str1);
                if tf==0 % Map doesn't have it
                    map(str1) = MapCnt;
                    Names(MapCnt) = cellstr(str1);
                    MapCnt=MapCnt+1;
                end
            end
            if Ui(1)==0
                Ui(1) = map(str1);
            else
                Ui = [Ui map(str1)];
            end
            str1='';i=i+1;
        else
            str1=strcat(str1,str(i));
        end
    end
    tf = isKey(map,str1);
    if tf==0 % Map doesn't have it
        map(str1) = MapCnt;
        Names(MapCnt) = cellstr(str1);
        MapCnt=MapCnt+1;
    end
    if Ui(1)==0
        Ui(1) = map(str1);
    else
        Ui = [Ui map(str1)];
    end
    % Job of user - Hashmap
    str = char(JobTitle(k));
    if JobMap.size == 0
        JobMap(str) = JobMapCnt;
        JobMapCnt=JobMapCnt+1;
    else
        tf = isKey(JobMap,str);
        if tf==0 % Map doesn't have it
            JobMap(str) = JobMapCnt;
            JobMapCnt=JobMapCnt+1;
        end
    end
    if tf==0
        UJi = JobMapCnt-1;
    else
        UJi = JobMap(str);
    end
    % Make all skill vectors in Ui contract to their weighted mean
    meanS = zeros(3,1);
    totCount=0;
    for it=1:length(Ui)
        meanS = meanS + (Scount(Ui(it))+1)*S(:,Ui(it));
        totCount = totCount + Scount(Ui(it))+1;
    end
    meanS = meanS/totCount;
    % Actually updating the skill vectors
    for it=1:length(Ui)
        S(:,Ui(it)) = ( (totCount - Scount(Ui(it)) - 1) * S(:,Ui(it)) + (Scount(Ui(it))+1)*meanS)/totCount;		%Cross weights taken into account
    end
    %Build the skill matrix Si for user i
    Si = zeros(3,length(Ui));
    for it=1:length(Ui)
        Si(:,it) = S(:,Ui(it));
        Scount(Ui(it))=Scount(Ui(it))+1;			%Also updating the Skill count
    end
    Bi = A*Si;          %Computing Bi, the matrix that is approximation to possible jobs
    if MarkJ(UJi)==0	% Job is new to the program
        MarkJ(UJi)=1;
        J(:,UJi) = mean(Bi')';
        Jcount(UJi) = 1;
    else				%Job is known to the program
        % Find nearest vector in Bi that matches J(Uji)
        min = Inf;
        minI = 1;
        for j=1:length(Bi)
            if min>norm(Bi(j)-J(UJi))
                min=norm(Bi(j)-J(UJi));
                minI = j;
            end
            % Adjust the job vector according to weighted mean
            J(:,UJi) = (Jcount(UJi)*J(:,UJi) + Bi(minI))/(Jcount(UJi)+1);
            Jcount(UJi)=Jcount(UJi)+1;
        end
    end
end
%% Plotting the Skill Space 
figure;
scatter3(S(1,1:MapCnt-1),S(2,1:MapCnt-1),S(3,1:MapCnt-1),'filled','b','DisplayName','Skill Vectors');    hold on;
legend('show');
%% Step 3 :- Given career goal of user i, find the related skills users will require
[ID2,T2] = xlsread('D1.xlsx');  % Read the query data
JobTitle2 = T2(2:length(T2),3);
Skills2 = T2(2:length(T2),2);

for k=1:length(ID2)
    str3=char(JobTitle2(k));  % Query that will be given
    disp('For ');disp(str3);
    Uji = JobMap(str3);       % Get the respective JobNo assigned
    b = J(:,Uji);
    x = inv(A)*b;             % A should be invertible
    min = Inf;
    minI = 1;
    for i=1:size(S, 2)
        if norm(x-S(:,i))<min
            min = norm(x-S(:,i));
            minI = i;
        end
    end
    for i=1:length(Names)
        if norm(S(:,minI)-S(:,i))<0.5*min     % Threshold is troublesome due to randomness
            disp(Names(i));
        end
    end
end