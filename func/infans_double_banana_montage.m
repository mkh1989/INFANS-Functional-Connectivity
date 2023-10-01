function newMontage = infans_double_banana_montage(oldMontage)
% this function converts referntial montage (Cz) to double banana montage (19 channels -> 18 channels)
% Input:
% 	oldMontage: a matrix containing eeg signals in referntial montage     (19 * samples)
%		ch1 : Fp1,  ch2 : Fp2,  ch3 : F3,  ch4 : F4,  ch5 : C3 
%		ch6 : C4,   ch7 : P3,   ch8 : P4,  ch9 : O1,  ch10: O2 
%		ch11: F7,   ch12: F8,   ch13: T3,  ch14: T4,  ch15: T5 
%		ch16: T6,   ch17: Fz,   ch18: Cz,  ch19: Pz
%
% Output:
%	newMontage: a matrix containing eeg signals in double banana montage (18 * samples)
%		ch1 : Fp2-F8, ch2 : F8-T4,   ch3 : T4-T6,  ch4 : T6-O2,  ch5 : O2-P4,  Ch6 : P4-C4
%		ch7 : C4-F4,  ch8 : F4-Fp2,  ch9 : Fz-Cz,  ch10: Cz-Pz,  ch11: Fp1-F7, Ch12: F7-T3
%		ch13: T3-T5,  ch14: T5-O1,   ch15: O1-P3,  ch16: P3-C3,  ch17: C3-F3,  Ch18: F3-Fp1

    newMontage = zeros(18,size(oldMontage,2));
    newMontage(1,:)  = oldMontage(2,:)  - oldMontage(12,:); % Fp2-F8
    newMontage(2,:)  = oldMontage(12,:) - oldMontage(14,:); % F8-T4
    newMontage(3,:)  = oldMontage(14,:) - oldMontage(16,:); % T4-T6
    newMontage(4,:)  = oldMontage(16,:) - oldMontage(10,:); % T6-O2
    newMontage(5,:)  = oldMontage(10,:) - oldMontage(8,:);  % O2-P4
    newMontage(6,:)  = oldMontage(8,:)  - oldMontage(6,:);  % P4-C4
    newMontage(7,:)  = oldMontage(6,:)  - oldMontage(4,:);  % C4-F4
    newMontage(8,:)  = oldMontage(4,:)  - oldMontage(2,:);  % F4-Fp2
    newMontage(9,:)  = oldMontage(17,:) - oldMontage(18,:); % Fz-Cz
    newMontage(10,:) = oldMontage(18,:) - oldMontage(19,:); % Cz-Pz
    newMontage(11,:) = oldMontage(1,:)  - oldMontage(11,:); % Fp1-F7
    newMontage(12,:) = oldMontage(11,:) - oldMontage(13,:); % F7-T3
    newMontage(13,:) = oldMontage(13,:) - oldMontage(15,:); % T3-T5
    newMontage(14,:) = oldMontage(15,:) - oldMontage(9,:);  % T5-O1
    newMontage(15,:) = oldMontage(9,:)  - oldMontage(7,:);  % O1-P3
    newMontage(16,:) = oldMontage(7,:)  - oldMontage(5,:);  % P3-C3
    newMontage(17,:) = oldMontage(5,:)  - oldMontage(3,:);  % C3-F3
    newMontage(18,:) = oldMontage(3,:)  - oldMontage(1,:);  % F3-Fp1

end