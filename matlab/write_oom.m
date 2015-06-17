function write_oom(filename, sig, tau, w0)
% write the given Oom in tom format to file
% tau must be anarray of size dim x dim x nO for an OOM or
% dim x dim x nO x nU for an IO-OOM
fid = fopen(filename, 'w');
nO = size(tau, 3);
if ndims(tau) == 3
    nU = 0;
else
    nU = size(tau, 4);
end
dim = size(tau,1);

fprintf(fid, '{"Type":"OOM","nU":%i,"nO":%i,"dim":%i,"sig":[[', nU, nO, dim);
for i = 1:dim-1
    fprintf(fid, '%0.15g,', sig(i));
end
fprintf(fid, '%0.15g]],"tau":', sig(dim));

fprintf(fid, '[')
for o = 1:nO
    fprintf(fid, '[')
    for u = 1:max(nU,1)
        fprintf(fid, '[');
        for i = 1:dim
            fprintf(fid, '[')
            for j = 1:dim
                if nU == 0
                    fprintf(fid, '%0.15g', tau(i,j,o));
                else
                    fprintf(fid, '%0.15g', tau(i,j,o,u));
                end
                if j == dim
                    fprintf(fid, ']')
                else
                    fprintf(fid, ',')
                end
            end
            if i == dim
                fprintf(fid, ']')
            else
                fprintf(fid, ',')
            end
        end
        if u == max(nU,1)
            fprintf(fid, ']')
        else
            fprintf(fid, ',')
        end
    end
    if o == nO
        fprintf(fid, ']')
    else
        fprintf(fid, ',')
    end
end

fprintf(fid, ',"w0":[');
for i = 1:dim-1
    fprintf(fid, '[%0.15g],', w0(i));
end
fprintf(fid, '[%0.15g]]}', w0(dim));
fclose(fid);
