clc
clear

%% basename of input file
basename = getenv('basename');
fprintf(basename);

%% path of matrix
path1 = getenv('prefix');
fprintf(path1);

%% the data type of matrix
helper = getenv('suffix');
fprintf(helper);

%% where to restore the results
path3 = getenv('storage');
fprintf(path3);

%% path of restore the gene density
path4 = getenv('pipe_path');
fprintf(path4);

%% chromosome numbers
chrom = 23;

%% bin number of smoothing (4 or 6)
smooth = 6;

%% resolutions
bin_size = getenv('bin_size');
fprintf(bin_size);
resolution = str2num(bin_size);

%% **********************below are calculting Compartment!!***************************
%% for each chromosome
genes_density = zeros(50,4);
ids1 =1;
ABcompartment = [];
gene = load('/bar/yliang/softwares/compartment_call_Bing_pipeline_2015/hg38_coding_gene_10kb.density');%this file contains the tss number in each 10kb bin on each chromosome
contributions = zeros(50,100);
split = load('/bar/yliang/softwares/compartment_call_Bing_pipeline_2015/chrom_split.txt');
for i = 1:chrom
    fprintf('chr%d is going to be processed!\n',i);
    %% obtain the interactions matrix
    %hh = [num2str(i-1),' ',num2str(i-1)];
    hh = num2str(i);
    data = load(strcat(path1,hh,helper));
    %% normalization for sequencing depth
    data(isnan(data)==1) = 0;
    data = data/sum(sum(data))*2000000;
    rows = size(data,1);
    %% from now on, we split the chromosome into pq arms
    if split(split(:,1)==i,2)~=0
        mid = split(i,2)/resolution;
        pos_part = [0,mid,rows];
        turn = 2;
        %adding begain, to split q arm into 2
        if i <=2
            mid2 = ceil(rows/2+mid/2);
            pos_part = [0,mid,mid2,rows];
            turn = 3;
        end
        %adding end
    else
        turn = 1;
        pos_part = [0,rows];
    end
    info_compart = [];
    for part = 1:turn
        %% zeros cols or rows positions
        datap = data(pos_part(part)+1:pos_part(part+1),pos_part(part)+1:pos_part(part+1));
        rowsp = size(datap,1);
        m_0 = sum(datap)==0;
        matrix_0 = ones(rowsp,rowsp);
        matrix_0(m_0==1,:)=0;
        matrix_0(:,m_0==1)=0;    
        %% smoothing data
        smooth_bin = ceil(smooth/2);
        sm_data1 = zeros(rowsp,rowsp);
        for j = 1:rowsp
            if m_0(j)~=1
                start1 = max(1,j-smooth_bin);
                over1 = min(rowsp,j+smooth_bin);
                filter = datap(:,start1:over1);
                filter(:,sum(filter,1)==0) = [];
                sm_data1(:,j) = mean(filter,2);
            end
        end
        sm_data2 = zeros(rowsp,rowsp);
        for j = 1:rowsp
            if m_0(j)~=1
                start1 = max(1,j-smooth_bin);
                over1 = min(rowsp,j+smooth_bin);
                filter = sm_data1(start1:over1,:);
                filter(sum(filter,2)==0,:) = [];
                sm_data2(j,:) = mean(filter,1);
            end
        end

        %% data normalization by distance
        expect_v = zeros(rowsp,1);
        for j = 1:rowsp
            diags = diag(datap,j);
            diags_h = diag(matrix_0,j);
            diags_u = diags(diags_h==1);
            expect_v(j) = mean(diags_u);
        end

        %% generate expect matrix
        expect = zeros(rowsp,rowsp);
        for j = 1:rowsp
            for k = 1:rowsp
                if m_0(j)~=1 && m_0(k)~=1 && j~=k
                    expect(j,k) = expect_v(abs(j-k));
                end
            end
        end
        %% smoothing expect
        sm_expect1 = zeros(rowsp,rowsp);
        for j = 1:rowsp
            if m_0(j)~=1
                start1 = max(1,j-smooth_bin);
                over1 = min(rowsp,j+smooth_bin);
                filter = expect(:,start1:over1);
                filter(:,sum(filter,1)==0) = [];
                sm_expect1(:,j) = mean(filter,2);
            end
        end
        sm_expect2 = zeros(rowsp,rowsp);
        for j = 1:rowsp
            if m_0(j)~=1
                start1 = max(1,j-smooth_bin);
                over1 = min(rowsp,j+smooth_bin);
                filter = sm_expect1(start1:over1,:);
                filter(sum(filter,2)==0,:) = [];
                sm_expect2(j,:) = mean(filter,1);
            end
        end
        %% observe/expect
        o_e = sm_data2./sm_expect2;
        o_e(isinf(o_e)==1) = 0;
        o_e(isnan(o_e)==1) = 0;
        %% add remove zeros
        p_0 = sum(o_e)~=0;
        o_e1 = o_e;%% restore the data
        o_e(sum(o_e)==0,:) = [];
        o_e(:,sum(o_e)==0) = [];
        %% correlation matrix
        corrs = corr(o_e,'type','Spearman');
        corrs(isinf(corrs)==1) = 0;
        corrs(isnan(corrs)==1) = 0;
        %% compute principle vectors and contributions
        [coeff,score,latent] = pca(corrs);
        result = coeff(:,1);
        contri = latent(1:100)/sum(latent);
        contributions(ids1,:) = contri';
        
        %% deal with the problems with zeros rows and cols
        result1 = result;
        result = zeros(rowsp,1);
        result(p_0) = result1;
        result = result*sqrt(sum(p_0));
        %% identifying the sign of vectors 
        gene_helper = zeros(rowsp,1);
        for j =pos_part(part)+1:pos_part(part+1)
            gene_helper(j-pos_part(part)) = sum(gene((j-1)*(resolution/10000)+1:j*(resolution/10000),i));
        end
        help1 = (result>0);
        help2 = (result<0);
        gene1 = sum(help1.*gene_helper)/sum(help1);
        gene2 = sum(help2.*gene_helper)/sum(help2);
        genes_density(ids1,:) = [i,gene1,gene2,max(gene1,gene2)/(gene1+gene2)];
        ids1 = ids1+1;
        if gene1<gene2
            result = -1*result;
        end
        result(p_0==0) = nan;
        info_compart = [info_compart',result']';
    end
    for j = 1:rows
        ss = resolution*(j-1)+1;
        oo = resolution*j;
        ABcompartment = [ABcompartment',[i,ss,oo,info_compart(j)]']';
    end

end
%% restore contributions
contributions = contributions(1:ids1-1,:);
dlmwrite(strcat(path3,'/', basename, '.contributions_each_chr_for_first100.txt'),contributions,...
    'delimiter','\t','precision',5);
dlmwrite(strcat(path3,'/', basename, '.AB_compartment_4_0_single_arm.txt'),ABcompartment,...
    'delimiter','\t','precision',9);
genes_density = genes_density(1:ids1-1,:);
dlmwrite(strcat(path3,'/', basename, '.gene_numbers_each_chrom_A_B.txt'),genes_density,...
    'delimiter','\t','precision',5);
%% finish !!!!