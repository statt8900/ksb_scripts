from db import sqlexecute,insertObject

commands = []
# commands.append("CREATE TABLE adsorbate (id integer primary key, created numeric, lastmodified numeric, name varchar, templatepath varchar)")
# commands.append("CREATE TABLE project (id integer primary key, created numeric, lastmodified numeric, name varchar)")
# commands.append("CREATE TABLE catalyst (id integer primary key, created numeric, lastmodified numeric, name varchar, metal varchar, crystalstructure varchar, facet varchar)")
# commands.append("CREATE TABLE catalystproperty (id integer primary key, created numeric, lastmodified numeric, catalystid integer, key varchar, value varchar, foreign key(catalystid) references catalyst(id))")
# commands.append("CREATE TABLE relaxation (id integer primary key, created numeric, lastmodified numeric, projectid integer, adsorbateid integer, catalystid integer, fullpath varchar, energy numeric, foreign key(adsorbateid) references adsorbate(id), foreign key(catalystid) references catalyst(id), foreign key(projectid) references project(id))")
# commands.append("CREATE TABLE configuration (id integer primary key, created numeric, lastmodified numeric, relaxationid integer, state varchar, aseadsorbateindex integer, asecatalystindex integer, foreign key(relaxationid) references relaxation(id))")
# commands.append("CREATE TABLE relaxationproperty (id integer primary key, created numeric, lastmodified numeric, relaxationid integer, key varchar, value varchar, foreign key(relaxationid) references relaxation(id))")
# commands.append("CREATE TABLE bindingenergy (id integer primary key, created numeric, lastmodified numeric, adsorbateid integer, catalystid integer, systemrelaxationid integer, cleanrelaxationid integer, electronicenergy numeric, freeenergy numeric, foreign key(cleanrelaxationid) references relaxation(id), foreign key(systemrelaxationid) references relaxation(id), foreign key(adsorbateid) references adsorbate(id), foreign key(catalystid) references catalyst(id))")
# commands.append("CREATE TABLE fbl (id integer primary key, created numeric, lastmodified numeric, relaxationid integer, fullpath varchar, tsenergy numeric, barrier numeric, foreign key(relaxationid) references relaxation(id))")
# commands.append("CREATE TABLE atom (id integer primary key, created numeric, atomicnumber integer, symbol varchar")
# commands.append("CREATE TABLE adsorbateatom (id integer primary key, created numeric, atomid integer, adsorbateid integer, foreign key(atomid) references atom(id), foreign key(adsorbateid) references adsorbate(id))")
# commands.append("CREATE TABLE catalystatom (id integer primary key, created numeric, atomid integer, catalystid integer, foreign key(atomid) references atom(id), foreign key(catalystid) references catalyst(id))")


commands.append("""CREATE TABLE bulk       
					(id        integer primary key
					,pth       varchar
					,name      varchar
					,stoich    varchar
					,bravais   varchar
					,inita     numeric
					,initb     numeric
					,initc     numeric
					,initalpha numeric
					,initbeta  numeric
					,initgamma numeric
					,initpos   varchar
					,mag       varchar
					,emt       numeric
					,comments  text)""") 

commands.append("""CREATE TABLE xc         
					(id     integer primary key
					,name   varchar
					,coeffs varchar)""") #how to store 8x8 matrix?

commands.append("""CREATE TABLE psp        
					(id    integer primary key
					,name  varchar
					,pth   varchar
					,nelec varchar)""")  #list where index = atomic number? can't print dict due to quotes

commands.append("""CREATE TABLE calc       
					(id       integer primary key
					,xc       varchar
					,psp      varchar
					,pw       integer
					,kpt      varchar
					,econv    numeric
					,fmax     numeric
					,xtol     numeric
					,mixing   numeric
					,nmix     integer
					,maxsteps integer
					,nbands   integer
					,sigma    numeric
					,magmom   integer)""")

commands.append("""CREATE TABLE bulkjob    
					(id           integer primary key
					,created      numeric
					,createdby    varchar
					,lastmodified numeric
					,kind         varchar
					,comments     text
					,status       varchar
					,bulk         varchar
					,calcid       integer
					,precalcxc    varchar
					,dft          varchar
					,timelim      numeric
					,foreign key(calcid) references calc(id)
					,foreign key(bulk) references bulk(pth))""")

commands.append("""CREATE TABLE bulkresult 
					(id          integer primary key
					,jobid       integer
					,a           numeric
					,b           numeric
					,c           numeric
					,alpha       numeric
					,beta        numeric
					,gamma       numeric
					,pos         varchar
					,magmom      varchar
					,energy      numeric
					,bulkmodulus numeric
					,bfit        numeric
					,xccoeffs    varchar
					,time        numeric
					,niter       integer
					,foreign key(jobid) references bulkjob(id))""")

for i,command in enumerate(commands): sqlexecute(command)