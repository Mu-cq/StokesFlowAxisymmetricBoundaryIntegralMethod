%rescale a function to get a desired integral value

function R = dfGetNormal(dfADD,df,nx,dfX)

df = df+dfADD;
dfx = df.*nx;
R = dfx-dfX;