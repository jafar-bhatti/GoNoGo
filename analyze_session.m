function [d_prime, data] = analyze_session(raw_data)

%NO_GO Data
fa = sum(raw_data == 3); %Number of false alarm trials
cr = sum(raw_data == 2); %Number of correct reject trials
no_go = fa + cr; %Number of no-go trials 

%GO data
hit = sum(raw_data == 1); %Number of hit trials
miss = sum(raw_data == 0); %Number of miss trials
go = hit + miss; %Number of go trials 

%Calculate hit rate 
hit_rate = hit/go;
        if hit_rate == 1 
            hit_rate = (hit-1)/go;
        end
        if hit_rate == 0
            hit_rate = 1/go;
        end
        
%Calculate fa rate 
fa_rate = fa/no_go;
        if fa_rate == 1
            fa_rate = (fa-1)/no_go;
        end
        if fa_rate == 0
            fa_rate = 1/no_go;

%Calculate d'
d_prime = norminv(hit_rate)-norminv(fa_rate);

%Summarize Data in a table 
Total_Trials = length(raw_data); %Total number of trials in the session
No_Go_Trials = no_go;
Go_Trials = go;
Hit_Trials = hit;
Miss_Trials = miss;
False_Alarm_Trials = fa;
Correct_Reject_Trials = cr;

data = table(Total_Trials, No_Go_Trials, Go_Trials, Hit_Trials, ... 
Miss_Trials, False_Alarm_Trials, Correct_Reject_Trials, d_prime);

end