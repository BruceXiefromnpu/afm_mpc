
function contents= freadf(fname)
fid = fopen(fname);

contents = '';
        while 1
            tline = fgets(fid);
            if ~ischar(tline), 
                break
            end
            contents = sprintf('%s%s', contents, tline);
            
        end
        
fclose(fid);
end


