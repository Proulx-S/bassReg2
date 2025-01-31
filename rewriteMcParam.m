function rewriteMcParam(prefix)
prefix = cellstr(prefix);
global srcAfni

ext = '.param.1D';

cmd = {srcAfni};
for i = 1:numel(prefix)
   [prefix{i} ext]
    cmd{end+1} = ['rm -f ' finalPreprocFiles.fTransCatList{r}];
    cmd{end+1} = ['head -1 ' strjoin(finalPreprocFiles.fTransList(r,1),' ') ' > ' finalPreprocFiles.fTransCatList{r}];
    cmd{end+1} = ['cat_matvec ' strjoin(flip(finalPreprocFiles.fTransList(r,:)),' ') ' >> ' finalPreprocFiles.fTransCatList{r}];
end

cellstr(prefix)
prefix

