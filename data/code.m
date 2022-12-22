clc
clear;
path = 'G:\Thesis\Modis\cj\C_re1\';
saveddir = 'G:\Thesis\Modis\cj\C_Result\'; % 图像新保存到的路径
%savedname = fullfile(saveddir,images(j).name); % 图像名称不变
subdir = dir(path);
for b=3:length(subdir)%第一层循环遍历子文件夹
    if(isequal(subdir(b),'.') || isequal(subdir(b),'..') || ~subdir(b).isdir)
        continue;
    end
    subdirpath = fullfile(path,subdir(b).name,'*.tif');
    images = dir(subdirpath); % 所有后缀为.tif的文件
    k=1;
    cd=2021-2012+1;
    L=length(images);
    new1=subdir(b).name;
    for c=1:length(images)%第二层循环遍历每个子文件夹中的图片
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
    geotiffwrite(savedname,result,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag)%注意修改路径
    disp('Sen slope finish!')
    
    for i=1:size(datasum,1)
        data=datasum(i,:);
        if min(data)>0       % 有效格点判定，我这里有效值在0以上
            sgnsum=[];
            for k=2:cd       %作用类似于sgn函数    xj-xi>0,sgn=1; xj-xi=0,sgn=0; xj-xi<0,sgn=-1;   (后减前)
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
                    sgnsum=[sgnsum;sgn];  %在sgnsum后面再加上sgn
                end
            end
            add=sum(sgnsum);
            sresult(i)=add;  %检验统计量S
        end
    end
    vars=cd*(cd-1)*(2*cd+5)/18;
    zc=zeros(m,n);
    sy=find(sresult==0);    %|Z|>1.96变化显著，|Z|<=1.96时变化不显著
    zc(sy)=0;                             %S=0时
    sy=find(sresult>0);
    zc(sy)=(sresult(sy)-1)./sqrt(vars);   %S>0时
    sy=find(sresult<0);
    zc(sy)=(sresult(sy)+1)./sqrt(vars);   %S<0时
    name1=strcat(new1,'_MKresult_C.tif');
    savedname1 = fullfile(saveddir,name1);
    geotiffwrite(savedname1,sresult,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag)%注意修改路径
    disp('MK finish!')
    
    
    result1=reshape(result,m*n,1);
    zc1=reshape(zc,m*n,1);
    tread =zeros(m,n);
    for i=1:size(datasum,1)
        
        if result1(i)>0   % Sen趋势B>0  上升
            if abs(zc1(i))>=2.58    %极显著上升
                tread(i)=4;
            elseif (1.96<=abs(zc1(i)))&&(abs(zc1(i))<2.58)      %显著上升
                tread(i)=3;
            elseif (1.645<=abs(zc1(i)))&&(abs(zc1(i))<1.96)     %微显著上升
                tread(i)=2;
            else		%不显著上升
                tread(i)=1;
            end
        elseif result1(i)<0  % Sen趋势B<0  下降
            if abs(zc1(i))>=2.58   %极显著下降
                tread(i)=-4;
            elseif (1.96<=abs(zc1(i)))&&(abs(zc1(i))<2.58)  %显著下降
                tread(i)=-3;
            elseif (1.645<=abs(zc1(i)))&&(abs(zc1(i))<1.96)   %微显著下降
                tread(i)=-2;
            else			%不显著下降
                tread(i)=-1;
            end
        else
            tread(i)=0;      % 无变化
        end
    end
    name2=strcat(new1,'_C显著性.tif');
    savedname2 = fullfile(saveddir,name2);
    geotiffwrite(savedname2,tread,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag)%注意修改路径
    disp('finish!')
    
end
disp('我是天才！！！！')