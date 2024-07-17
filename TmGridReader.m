function TmGrid=TmGridReader(path)
pattern = ["h00","h06","h12","h18"];
TmGrid=[];
TF = endsWith(path,pattern);
if(find(TF==1)==0)
    return;
end
fid=fopen(path,'rt');
if (fid==-1)
    %msgbox('file invalide','warning','warn');
    return;
end

temp=[];
while ~feof(fid)
    s=fgets(fid);
    sizes=strlength(s);
    if(sizes<2)
        continue;
    end
    str=strtrim(s);
    if(strcmp(str(1:5),'90.00'))
        continue;
    end
    str=strsplit(str,' ');

    if(sizes>10)
        temp=[temp str2double(str)];
    end
    if(sizes<10)
        TmGrid=[TmGrid temp];
        temp=[];
        continue;
    end
end
fclose(fid);
