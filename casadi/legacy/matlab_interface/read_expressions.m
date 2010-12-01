function ret=read_expressions(path_to_file)
resfile = fopen(path_to_file);

ret = struct;
while 1
    % One variable per line
    str = fgetl(resfile);
    if ~ischar(str), break, end

    [name,count,errmsg,nextindex] = sscanf(str,'%s',1);
    str = str(nextindex:end);
    
    % Try to read a matrix
    [d,count,errmsg,nextindex] = sscanf(str,' [%f,%f]',2);
 
 if(length(errmsg) == 0)
    str = str(nextindex:end);

    form = repmat('%f,',1,d(2));
    form(end) = ')';
    form = ['(',form];
    form = repmat([form,','],1,d(1));
    form(end) = ')';
    form = ['(',form];
 
    el = nan([d(2),d(1)]);
    el(:) = sscanf(str,form,d(1)*d(2));
   
    ret = setfield(ret,name,el');
    continue;    
 end
     % Try to read a vector
    [d,count,errmsg,nextindex] = sscanf(str,' [%f]',1);
     
    if(length(errmsg) == 0)
        str = str(nextindex:end);

        form = repmat('%f,',1,d);
        form(end) = ')';
        form = ['(',form];

        el = nan(d,1);
        el(:) = sscanf(str,form,d);
        
        ret = setfield(ret,name,el);
        
       continue;
    end
    
    % Read scalar
    [d,count,errmsg,nextindex] = sscanf(str,' %f',1);
    assert(length(errmsg) == 0);
    ret = setfield(ret,name,d);
end
fclose(resfile);
disp (sprintf('successfully read %s', path_to_file))
