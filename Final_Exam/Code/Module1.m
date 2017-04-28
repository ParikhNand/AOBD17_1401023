% AOBD FInal Exam Paper - Career Recommender System from LinkedIn's user's profiles 
% Module 1 -> Suggesting skills that are similar to the current skills
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

%% Step 2 :- Given skills of user i, find other related skills
[ID2,T2] = xlsread('D1.xlsx');  % Read the query data
JobTitle2 = T2(2:length(T2),3);
Skills2 = T2(2:length(T2),2);

for k=1:length(ID2) % For each query
    Ui = zeros(1,1);
    % Get the corresponding skill's value from HashMap
    str = char(Skills2(k));
    str2 = '';
    for i=1:length(str)
        if str(i)==','
            if Ui(1)==0
                Ui(1) = cell2mat(values(map,cellstr(str2)));
            else
                val = cell2mat(values(map,cellstr(str2)));
                Ui = [Ui val];
            end
            str2='';i=i+1;
        else
            str2=strcat(str2,str(i));
        end
    end
    if Ui(1)==0
        Ui(1) = cell2mat(values(map,cellstr(str2)));
    else
        val = cell2mat(values(map,cellstr(str2)));
        Ui = [Ui val];
    end
    % Create a Skill matrix X
    X = Ui;
    Xnew = zeros(1, 2*size(X,2));	%For each vector in X, we find 2 nearest skill vectors
    for i=1:size(X,2)
        %find nearest 2 skills
        min1 = inf;
        min1skill = 1;
        min2 = inf;
        min2skill = 2;
        for j=1:length(Names)
            if abs(norm(S(:,X(:,i))-S(:,j)))<min1 && j~=X(:,i)
                min2 = min1;
                min2skill = min1skill;
                min1 = abs(norm(S(:,X(:,i))-S(:,j)));
                min1skill = j;
            else if abs(norm(S(:,X(:,i))-S(:,j)))<min2  && j~=X(:,i)
                    min2 = abs(norm(S(:,X(:,i))-S(:,j)));
                    min2skill = j;
                end
            end
            Xnew(1,2*i-1) = min1skill;
            Xnew(1,2*i) =  min2skill;
        end
    end
    %Sorting and displaying "unique" results
    sort(Xnew);
    disp(['Suggested skills for user ',num2str(k)]);
    disp(Names(Xnew(1))); disp(Names(Xnew(2)));
    S_sug = zeros(3,500); cnt=2;
    S_sug(:,1) = S(:,Xnew(1));
    for i=2:size(Xnew, 2)
        if Xnew(i)~=Xnew(i-1)
            S_sug(:,cnt) = S(:,Xnew(i));    cnt=cnt+1;
        end
    end
end
% Mark the suggested skills differently and show in figure
scatter3(S_sug(1,1:MapCnt-1),S_sug(2,1:MapCnt-1),S_sug(3,1:MapCnt-1),'*','r','DisplayName','Suggested Skills');
legend('show'); title('Skill Space');