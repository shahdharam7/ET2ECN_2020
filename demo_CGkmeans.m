%% demo_CG-kmeans

clc;
close all;
clear all;

%% Image Read
s=load('Cuprite.mat');  % link for data source : https://rslab.ut.ac.ir/data
p=s.nRow;
q=s.nCol;
Bands=188;
Y=s.Y;
x=hyperConvert3d(Y,p,q,Bands);
for i=1:1:Bands
    R=x(:,:,i);
    minr=min(min(R));
    maxr=max(max(R));
    R1=(R-minr)/(maxr-minr);
    R1(R1==0)=0.0001;
    xp(:,:,i)=R1;
end

%% Virtual Dimension
VD=12;

%% CG-kmeans algorithm
[endmemberindex] = CGkmeans(Y,VD);
endmemberindex_CGKMEANS=change_index(endmemberindex,p,q);

%% VCA algorithm
[U_VCA,e_index,snrEstimate]=hyperVca(Y,VD);
endmemberindex_VCA=change_index(e_index,p,q);    

%% GT 
t1=load('groundTruth_Cuprite_nEnd12.mat');
gt=t1.M;
tit=t1.cood;
n1=gt(3:103,:);
n2=gt(114:147,:);
n3=gt(168:220,:);
gt=[n1;n2;n3];
[gt_m,gt_n]=size(gt);

%% Total Spectral Angle Mapper (TSAM) calculations
for i=1:gt_n
    for j=1:Bands
        extracted_CGKMEANS(j,i)=xp(endmemberindex_CGKMEANS(i,1),endmemberindex_CGKMEANS(i,2),j);
        extracted_VCA(j,i)=xp(endmemberindex_VCA(i,1),endmemberindex_VCA(i,2),j);
    end
end

[ex_m,ex_n]=size(extracted_VCA);
store_CGKMEANS=[0,0];
store_VCA=[0,0];
sam_VCA=0;
sam_CGKMEANS=0;
sam_total_CGKMEANS=0;
sam_total_VCA=0;

for i=1:gt_n
    for j=1:ex_n
            Mat_SAM_VCA(i,j)=real(acos(dot(gt(:,i),extracted_VCA(:,j))/(norm(gt(:,i)*norm(extracted_VCA(:,j))))));
            Mat_SAM_CGKMEANS(i,j)=real(acos(dot(gt(:,i),extracted_CGKMEANS(:,j))/(norm(gt(:,i)*norm(extracted_CGKMEANS(:,j))))));
    end
end

for i=1:gt_n
    %VCA
    [max_value1,mrow]=min(Mat_SAM_VCA);
    [max_value,col_VCA]=min(max_value1);
    sam_total_VCA=sam_total_VCA+max_value;
    sam_VCA=[sam_VCA;max_value];
    row_VCA=mrow(col_VCA);
    s1=[row_VCA,col_VCA];
    store_VCA=[store_VCA;s1];
    save_VCA(row_VCA)=max_value;
    Mat_SAM_VCA(row_VCA,:)=[100*ones];
    Mat_SAM_VCA(:,col_VCA)=[100*ones];
    %CGKMEANS
    [max_value1,mrow]=min(Mat_SAM_CGKMEANS);
    [max_value,col_CGKMEANS]=min(max_value1);
    sam_total_CGKMEANS=sam_total_CGKMEANS+max_value;
    sam_CGKMEANS=[sam_CGKMEANS;max_value];
    row_CGKMEANS=mrow(col_CGKMEANS);
    s1=[row_CGKMEANS,col_CGKMEANS];
    store_CGKMEANS=[store_CGKMEANS;s1];
    save_CGKMEANS(row_CGKMEANS)=max_value;
    Mat_SAM_CGKMEANS(row_CGKMEANS,:)=[100*ones];
    Mat_SAM_CGKMEANS(:,col_CGKMEANS)=[100*ones];
end

rms_sae=[rms(save_CGKMEANS);
    rms(save_VCA)];
rms_sae = radtodeg(rms_sae);

disp('RMSSAE of VCA');
disp(rms_sae(2));
disp('RMSSAE of CGKMEANS');
disp(rms_sae(1));