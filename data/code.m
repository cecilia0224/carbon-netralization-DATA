clc
clear;
path = 'G:\Thesis\Modis\cj\C_re1\';
saveddir = 'G:\Thesis\Modis\cj\C_Result\'; % ͼ���±��浽��·��
%savedname = fullfile(saveddir,images(j).name); % ͼ�����Ʋ���
subdir = dir(path);
for b=3:length(subdir)%��һ��ѭ���������ļ���
    if(isequal(subdir(b),'.') || isequal(subdir(b),'..') || ~subdir(b).isdir)
        continue;
    end
    subdirpath = fullfile(path,subdir(b).name,'*.tif');
    images = dir(subdirpath); % ���к�׺Ϊ.tif���ļ�
    k=1;
    cd=2021-2012+1;
    L=length(images);
    new1=subdir(b).name;
    for c=1:length(images)%�ڶ���ѭ������ÿ�����ļ����е�ͼƬ
        filename1 = fullfile(path,subdir(b).name,images(1).name);
        filename = fullfile(path,subdir(b).name,images(c).name);
        disp(filename);
        [a,R]=geotiffread(filename1);
        info=geotiffinfo(filename1);
        [m,n]=size(a);
        %cd=2021-2012+1
        datasum=zeros(m*n,11)+NaN;
    end
    for d=1:L
        file = fullfile(path,subdir(b).name,images(d).name);
        data=importdata(file);
        disp('yeah!')
        
        data=reshape(data,m*n,1);
        
        if (k>11)
            k=1;
        else
            datasum(:,k)=data;
            k=k+1;
            disp('+1')
        end
    end
    result=zeros(m,n)+NaN;
    sresult=zeros(m,n)+NaN;
    
    for i=1:size(datasum,1)
        data=datasum(i,:);
        if min(data)>0
            valuesum=[];
            for k1=2:cd
                for k2=1:(k1-1)
                    %                         disp(data(k1));
                    %                         disp(data(k2));
                    cz=data(k1)-data(k2);
                    %                         disp(cz);
                    j=k1-k2;
                    value=cz./j;
                    valuesum=[valuesum;value];
                end
            end
            value=median(valuesum);
            result(i)=value;
        end
        %
    end
    
    name=strcat(new1,'_trend_C.tif');
    savedname = fullfile(saveddir,name);
    geotiffwrite(savedname,result,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag)%ע���޸�·��
    disp('Sen slope finish!')
    
    for i=1:size(datasum,1)
        data=datasum(i,:);
        if min(data)>0       % ��Ч����ж�����������Чֵ��0����
            sgnsum=[];
            for k=2:cd       %����������sgn����    xj-xi>0,sgn=1; xj-xi=0,sgn=0; xj-xi<0,sgn=-1;   (���ǰ)
                for j=1:(k-1)
                    sgn=data(k)-data(j);
                    if sgn>0
                        sgn=1;
                    else
                        if sgn<0
                            sgn=-1;
                        else
                            sgn=0;
                        end
                    end
                    sgnsum=[sgnsum;sgn];  %��sgnsum�����ټ���sgn
                end
            end
            add=sum(sgnsum);
            sresult(i)=add;  %����ͳ����S
        end
    end
    vars=cd*(cd-1)*(2*cd+5)/18;
    zc=zeros(m,n);
    sy=find(sresult==0);    %|Z|>1.96�仯������|Z|<=1.96ʱ�仯������
    zc(sy)=0;                             %S=0ʱ
    sy=find(sresult>0);
    zc(sy)=(sresult(sy)-1)./sqrt(vars);   %S>0ʱ
    sy=find(sresult<0);
    zc(sy)=(sresult(sy)+1)./sqrt(vars);   %S<0ʱ
    name1=strcat(new1,'_MKresult_C.tif');
    savedname1 = fullfile(saveddir,name1);
    geotiffwrite(savedname1,sresult,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag)%ע���޸�·��
    disp('MK finish!')
    
    
    result1=reshape(result,m*n,1);
    zc1=reshape(zc,m*n,1);
    tread =zeros(m,n);
    for i=1:size(datasum,1)
        
        if result1(i)>0   % Sen����B>0  ����
            if abs(zc1(i))>=2.58    %����������
                tread(i)=4;
            elseif (1.96<=abs(zc1(i)))&&(abs(zc1(i))<2.58)      %��������
                tread(i)=3;
            elseif (1.645<=abs(zc1(i)))&&(abs(zc1(i))<1.96)     %΢��������
                tread(i)=2;
            else		%����������
                tread(i)=1;
            end
        elseif result1(i)<0  % Sen����B<0  �½�
            if abs(zc1(i))>=2.58   %�������½�
                tread(i)=-4;
            elseif (1.96<=abs(zc1(i)))&&(abs(zc1(i))<2.58)  %�����½�
                tread(i)=-3;
            elseif (1.645<=abs(zc1(i)))&&(abs(zc1(i))<1.96)   %΢�����½�
                tread(i)=-2;
            else			%�������½�
                tread(i)=-1;
            end
        else
            tread(i)=0;      % �ޱ仯
        end
    end
    name2=strcat(new1,'_C������.tif');
    savedname2 = fullfile(saveddir,name2);
    geotiffwrite(savedname2,tread,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag)%ע���޸�·��
    disp('finish!')
    
end
disp('������ţ�������')