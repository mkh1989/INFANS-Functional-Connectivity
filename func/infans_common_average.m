function newMontage = infans_common_average(oldMontage)
% this function converts the referntial montage to the common average Montage
newMontage = double(oldMontage - mean(oldMontage));
end