function [ F1score, recall, precision, TP, FP, FN,missed_saccades,onset_lag,offset_lag,extra_saccades] =...
    EventConfusion_Expanded( reference, predicted )
% [ F1score, recall, precision, TP, FP, FN ] = EventConfusion( reference, predicted )
% reference: Sleep-spindle detection: crowdsourcing and evaluating
% performance of experts, non-experts and automated methods 
% Simon C Warby, etc. 2014
%
% Code further modified From 'A parametric model for saccadic eye movement'
% IEEE x_signal Processing in Medicine and Biology Symposium (SPMB), December 2016.
% Modified by Seth Konig 3/14/17, added onset and offset_lag
%Modified by Ben Voloh 4/17/18, intersection/union was slow, optimized

%BV: 04/17/2018: check theyre both teh same size for later things to work
if isrow(reference); reference = reference'; end
if isrow(predicted); predicted = predicted'; end

Threshold = 0.20; 
[E, I] = bwlabel(reference);
[D, J] = bwlabel(predicted);
O_ED = zeros(I,J);
TP_candidate = zeros(I,J);
Event_match = zeros(I,J);
Detection_match = zeros(I,J);
TP = zeros(I,J);
onset_lag = NaN(1,I);
offset_lag = NaN(1,I);

FN = zeros(I,1);
FP = zeros(J,1);

for i = 1:I
    for j = 1:J
        temp_intersect = E==i & D==j; %BV: 04/17/2018, intersect/union is slow
        if any(temp_intersect) 
            temp_union = E==i | D==j;
            O_ED(i,j) = sum(temp_intersect)/sum(temp_union);
        end
    end
end

for i = 1:I
    for j = 1:J
        if O_ED(i,j) > Threshold
            TP_candidate(i,j) = 1;
        else
            TP_candidate(i,j) = 0;
        end
    end
end

for i = 1:I
    if sum(TP_candidate(i,:)) > 0
        [~, idx] = max(O_ED(i,:));
        Event_match(i,idx) = 1;
    elseif sum(O_ED(i,:)) == 0
        FN(i) = 1;
    end
end

for j = 1:J
    if sum(TP_candidate(:,j)) > 0
        [~, idx] = max(O_ED(:,j));
        Detection_match(idx,j) = 1;
    elseif sum(O_ED(:,j)) == 0
        FP(j) = 1;  
        %add FP index so can calculate duration and amplitude of false positive
    end
end
Best_match = Event_match + Detection_match;

for i = 1:I
    for j = 1:J
        if Best_match(i,j) == 2
            TP(i,j) = 1;
            Best_match(i,:) = 0;
            Best_match(:,j) = 0;
            
            true_time_pts = find(E==i);
            predicted_time_pts = find(D==j);
            onset_lag(i) = true_time_pts(1)-predicted_time_pts(1);
            offset_lag(i) = true_time_pts(end)-predicted_time_pts(end);
        end
    end
end

if sum(sum(Best_match)) ~= 0
    O_ED2 = zeros(I,J);
    O_ED2(Best_match == 1) = 1;
    TP_candidate2 = O_ED2;
    Event_match2 = zeros(I,J);
    Detection_match2 = zeros(I,J);
    for i = 1:I
        if sum(TP_candidate2(i,:)) > 0
            idx = max(O_ED2(i,:));
            Event_match2(i,idx) = 1;
        end
    end
    for j = 1:J
        if sum(TP_candidate2(:,j)) > 0
            idx = max(O_ED2(:,j));
            Detection_match2(idx,j) = 1;
        end
    end
    Best_match2 = Event_match2 + Detection_match2;
    for i = 1:I
        for j = 1:J
            if Best_match2 == 2
                TP(i,j) = 1;
            end
        end
    end
end

for i = 1:I
    if sum(TP(i,:)) == 0
        FN(i) = 1;
    end
end

for j = 1:J
    if sum(TP(:,j)) == 0
        FP(j) = 1;
    end
end

missed_saccades = FN;
extra_saccades = FP;
TP = sum(sum(TP));
FN = sum(FN);
FP = sum(FP);

recall = TP/(TP+FN);  %aka sensitivity
precision = TP/(TP+FP); %aka positive predictive value
F1score = 2*(precision*recall)/(precision+recall); %measure of accuracy    