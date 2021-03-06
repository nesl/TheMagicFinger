run the following files with 2345
run_singlept_manual.m
run_singlept_manual_theta_out.m
angle_vs_acc_140cm.m
angle_vs_acc_75cm.m

with theta :
Type 2 : RF : Accuracy: 0.841463
linear Accuracy: 0.463415
rbf Accuracy: 0.664634
quadratic Accuracy: 0.567073

Type 1 : RF Accuracy: 0.810976
linear Accuracy: 0.451220
rbf Accuracy: 0.646341
quadratic Accuracy:0.500000, 0.57926

wihtout theta
Type 2 : RF Accuracy: 0.804878
linear Accuracy: 0.40853,0.432927
rbf Accuracy: 0.518293,0.554878
quadratic Accuracy: 0.500,0.469512, 0.481707

Type 1: RF Accuracy:0.780488, 0.798780, 0.829268,0.835366
linear Accuracy: 0.414634
rbf Accuracy:  0.506098
quadratic Accuracy: 0.500 ,0.487805
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
type 2 with theta
% ===== Non-cooperative detection =====
% Training on 375, Testing on 160
% ===== SVM linear=====
% Precision: 0.533333
% Recall: 0.923077
% F score: 0.676056
% ===== SVM rbf=====
% Precision: 0.619048
% Recall: 0.962963
% F score: 0.753623
% ===== SVM poly =====
% Precision: 0.666667
% Recall: 0.962963
% F score: 0.787879
% ===== RF =====
% Precision: 0.657895
% Recall: 0.925926
% F score: 0.769231
% ===== Min-cooperative detection =====
% linear Accuracy: 0.580247
%rbf Accuracy: 0.765432
%poly Accuracy: 80.2
% RF Accuracy: 0.93
%%%%%%%%%%%%%%%%%%%%%
type 1 rij with theta
% ===== SVM linear=====
% Precision: 0.418182
% Recall: 0.851852
% F score: 0.560976
% ===== SVM rbf=====
% Precision: 0.590909
% Recall: 0.928571
% F score: 0.722222
% ===== SVM poly=====
% Precision: 0.593750
% Recall: 0.826087
% F score: 0.690909
% ===== RF =====
% Precision: 0.471698
% Recall: 0.925926
% F score: 0.625000
% ===== RF Min-cooperative detection =====
% Accuracy: 0.913580
% linear Accuracy:  0.567901
% reb Accuracy: 0.740741
% 'poly' Accuracy: 0.728395
%%----------------------------------

%without theta type 2
% ===== Non-cooperative detection =====
% Training on 375, Testing on 160
% ===== SVM linear=====
% Precision: 0.433962
% Recall: 0.884615
% F score: 0.582278
% ===== SVM rbf=====
% Precision: 0.526316
% Recall: 0.833333
% F score: 0.645161
% ===== SVM poly=====
% Precision: 0.571429
% Recall: 0.869565
% F score: 0.689655
% 
% ===== RF =====
% Precision: 0.588235
% Recall: 0.869565
% F score: 0.701754
% =====RF Min-cooperative detection =====
% Accuracy: 0.913580
% ===== SVM linear==
% Accuracy: 0.555556
% ===== SVM rbf=====
% Accuracy: 0.703704
% ===== SVM poly==== 
% Accuracy: 0.740741


----------------------
wihtout theta rij

===== SVM linear=====
Precision: 0.415094
Recall: 0.846154
F score: 0.556962

====== SVM rbf ====
Precision: 0.434783
Recall: 0.869565
F score: 0.579710
===== SVM poly=====
Precision: 0.476190
Recall: 0.800000
F score: 0.597015
===== RF =====
Precision: 0.564103
Recall: 0.814815
F score: 0.666667

===== linear Min-cooperative detection =====
Accuracy: 0.567
rbf Accuracy: 0.691358
===== poly Min-cooperative detection =====
Accuracy: 0.716049
===== RF Min-cooperative detection =====
Accuracy: 0.864198
