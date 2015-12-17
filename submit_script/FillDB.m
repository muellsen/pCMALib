function [something] = FillDB(name,outputs,repeats,indbname,nrproc)
        
        %% DB Stuff
        % Database Server
        host = 'localhost';

        % Database Username/Password
        user = 'root';
        password = 'montecarlo';

        % Database Name
        %dbName = 'mysql'; 
        dbName = 'mysql'; 

        % JDBC Parameters
        jdbcString = sprintf('jdbc:mysql://%s/%s', host, dbName);
        jdbcDriver = 'com.mysql.jdbc.Driver';

        % Set this to the path to your MySQL Connector/J JAR
        javaaddpath('/Users/ofgeorg/Documents/mysql-connector-java-5.1.12/mysql-connector-java-5.1.12-bin.jar')
        %javaaddpath('C:\Users\Georg\Downloads\mysql-connector-java-5.1.10\mysql-connector-java-5.1.10\mysql-connector-java-5.1.10-bin.jar')
        path(path,'/Users/ofgeorg/Documents/Submit Scrips Brutus')
        % Create the database connection object
            dbConn = database(dbName, user , password, jdbcDriver, jdbcString);

            % Check to make sure that we successfully connected
            if isconnection(dbConn)
                % Fetch the symbol, market cap, and last close for the 10 largest
            % market cap ETFs
            % If the connection failed, print the error message
            %%here the important stuff starts

            %first check if db exists
            string = ['CREATE DATABASE IF NOT EXISTS ' indbname];
            res = exec(dbConn, string);
            close(dbConn); 
            else
                disp(sprintf('Connection failed: %s', dbConn.Message));
            end
        %then create the tables according to the first file
        dbName = indbname; 
        
        jdbcString = sprintf('jdbc:mysql://%s/%s', host, dbName);
        jdbcDriver = 'com.mysql.jdbc.Driver';
        dbConn = database(dbName, user , password, jdbcDriver, jdbcString);
        % Check to make sure that we successfully connected
        if isconnection(dbConn)
        
        %create settings table
        [varnames,varvalue,vartypes,count] = Readsettings('output1repeat1',nrproc); %read first file        
        string = AddLine('','CREATE TABLE IF NOT EXISTS `setting` (');         
        string = AddLine(string,'`id` int(11) NOT NULL AUTO_INCREMENT,');
        for i = 1:count           
            string = AddLine(string,['`' varnames{i} '` ' vartypes{i} ' NOT NULL,']);
        end
        string = AddLine(string,'PRIMARY KEY (`id`)');
        string = AddLine(string,') ENGINE=MyISAM  DEFAULT CHARSET=latin1;');
        res = exec(dbConn, string);
        cd ('..');
        [status, message, messageid] = rmdir('output1repeat1','s');
        
        %create result table
        
        string = 'CREATE TABLE IF NOT EXISTS `result` (';
        string = AddLine(string,'`id` int(11) NOT NULL AUTO_INCREMENT,');
        string = AddLine(string,'`setting_id` int(11) NOT NULL,');
        if (nrproc > 1)
            string = AddLine(string,'`proc` int(11) NOT NULL,');
        end
        string = AddLine(string,'`out_bestever_f` double NOT NULL,');
        string = AddLine(string,'`out_bestever_evals` double NOT NULL,');
        string = AddLine(string,'`out_countEval` double NOT NULL,');
        string = AddLine(string,'`out_countOutOfBounds` double NOT NULL,');
        string = AddLine(string,'`out_countEvalNaN` double NOT NULL,');
        string = AddLine(string,'`out_countIter` double NOT NULL,');
        string = AddLine(string,'`out_lambda` double NOT NULL,');
        string = AddLine(string,'`out_mueff` double NOT NULL,');
        string = AddLine(string,'`out_mu` double NOT NULL,');
        string = AddLine(string,'`out_N` double NOT NULL,');
        string = AddLine(string,'`seed` double NOT NULL,');
        string = AddLine(string,'`out_stopflag` varchar(40) NOT NULL,');
        string = AddLine(string,'PRIMARY KEY (`id`)');
        string = AddLine(string,') ENGINE=MyISAM  DEFAULT CHARSET=latin1;');
        res = exec(dbConn, string);
        
        %create x start table
        string = 'CREATE TABLE IF NOT EXISTS `xstart` (';
        string = AddLine(string,'`id` int(11) NOT NULL AUTO_INCREMENT,');
        string = AddLine(string,'`setting_id` int(11) NOT NULL,');
        string = AddLine(string,'`x` double NOT NULL,');        
        string = AddLine(string,'PRIMARY KEY (`id`)');
        string = AddLine(string,') ENGINE=MyISAM  DEFAULT CHARSET=latin1;');
        res = exec(dbConn, string);
        
        %create bestever_x start table
        string = 'CREATE TABLE IF NOT EXISTS `bestever_x` (';
        string = AddLine(string,'`id` int(11) NOT NULL AUTO_INCREMENT,');
        string = AddLine(string,'`setting_id` int(11) NOT NULL,');
        string = AddLine(string,'`x` double NOT NULL,');        
        string = AddLine(string,'PRIMARY KEY (`id`)');
        string = AddLine(string,') ENGINE=MyISAM  DEFAULT CHARSET=latin1;');
        res = exec(dbConn, string);
        
        %create weights table
        
        string = 'CREATE TABLE IF NOT EXISTS `weights` (';
        string = AddLine(string,'`id` int(11) NOT NULL AUTO_INCREMENT,');
        string = AddLine(string,'`setting_id` int(11) NOT NULL,');
        string = AddLine(string,'`weight` double NOT NULL,');        
        string = AddLine(string,'PRIMARY KEY (`id`)');
        string = AddLine(string,') ENGINE=MyISAM  DEFAULT CHARSET=latin1;');
        res = exec(dbConn, string);
        
        
        %now we start filling it
        for i = 1:outputs
            for j = 1:repeats
                for p = 1:nrproc
                if (nrproc < 2)
                    proc_string = '';
                else
                    proc_string = ['_' num2str(p-1)];
                end
                    
                %insert the complete setting
                [varnames,varvalue,vartypes,count] = Readsettings(['output' num2str(i) 'repeat' num2str(j)],nrproc);
                string = AddLine('',['INSERT INTO `' indbname '`.`setting` (`id`']);
                for k = 1:count
                    string = AddLine(string,[',`' varnames{k} '`']);
                end
                string = AddLine(string,[') VALUES (NULL']);
                for k = 1:count
                    string = AddLine(string,[', ''' varvalue{k} '''']);
                end
                string = AddLine(string,[');']);
                res = exec(dbConn, string);
                string = ['SELECT LAST_INSERT_ID();']; %get the auto id
                result = get(fetch(exec(dbConn, string )), 'Data');
                id = result{1};
                
                %insert the results
                out_bestever_evals = load(['out_bestever_evals' proc_string '.txt']);                             
                out_bestever_f = load([ 'out_bestever_f' proc_string '.txt']);
                out_bestever_x = load([ 'out_bestever_x' proc_string '.txt']);
                %load 'out_BFGSexitcounter.txt'
                out_countEvalNaN = load([ 'out_countEvalNaN' proc_string '.txt']);
                out_countEval = load([ 'out_countEval' proc_string '.txt']);
                out_countIter = load([ 'out_countIter' proc_string '.txt']);
                out_countOutOfBounds = load([ 'out_countOutOfBounds' proc_string '.txt']);
                %load 'out_funcName.txt'
                out_insigma = load([ 'out_insigma' proc_string '.txt']);
                out_lambda = load([ 'out_lambda' proc_string '.txt']);
                out_mueff = load([ 'out_mueff' proc_string '.txt']);
                out_mu = load([ 'out_mu' proc_string '.txt']);
                out_N = load([ 'out_N' proc_string '.txt']);
                %load 'out_settings.txt'
                fid = fopen(['out_stopflag' proc_string '.txt'],'r');
                out_stopflag = fscanf(fid,'%s' );
                fclose(fid);
                %load 'out_stopflag.txt'
                %load 'out_weights.txt'
                out_xstart = load([ 'out_xstart' proc_string '.txt']);
                seed = load([ 'seed' proc_string '.txt']);
                %string = ['INSERT INTO `' name '`.`result` (`id`,
                %`setting_id`, `f`, `evals`, `countEval`, `countOutOfBounds`,`countEvalNaN`,`out_countIter`,`out_lambda`,`out_mueff`,`out_mu`,`out_N`,`seed`,`stopflag`) VALUES (NULL, ''' num2str(id) ''', ''' num2str(out_bestever_f) ''', ''' num2str(out_bestever_evals) ''', ''' num2str(out_countEval) ''', ''' num2str(out_countOutOfBounds) ''' , ''' num2str(stopflag) ''');'];
                string = ['INSERT INTO `' indbname '`.`result` (`id`, `setting_id`'];
                if (nrproc > 1)
                    string = AddLine(string,',`proc` ');
                end
                string = AddLine(string,',`out_bestever_f` ');
                string = AddLine(string,',`out_bestever_evals` ');
                string = AddLine(string,',`out_countEval` ');
                string = AddLine(string,',`out_countOutOfBounds` ');
                string = AddLine(string,',`out_countEvalNaN` ');
                string = AddLine(string,',`out_countIter` ');
                string = AddLine(string,',`out_lambda` ');
                string = AddLine(string,',`out_mueff` ');
                string = AddLine(string,',`out_mu` ');
                string = AddLine(string,',`out_N` ');
                string = AddLine(string,',`seed` ');
                string = AddLine(string,',`out_stopflag` ');
                string = AddLine(string,') VALUES (NULL');
                string = AddLine(string,[',''' num2str(id) '''']);
                if (nrproc > 1)
                    string = AddLine(string,[',''' num2str(p-1) '''']);
                end
                string = AddLine(string,[',''' num2str(out_bestever_f) '''']);
                string = AddLine(string,[',''' num2str(out_bestever_evals) '''']);
                string = AddLine(string,[',''' num2str(out_countEval) '''']);
                string = AddLine(string,[',''' num2str(out_countOutOfBounds) '''']);
                string = AddLine(string,[',''' num2str(out_countEvalNaN) '''']);
                string = AddLine(string,[',''' num2str(out_countIter) '''']);
                string = AddLine(string,[',''' num2str(out_lambda) '''']);
                string = AddLine(string,[',''' num2str(out_mueff) '''']);
                string = AddLine(string,[',''' num2str(out_mu) '''']);
                string = AddLine(string,[',''' num2str(out_N) '''']);
                string = AddLine(string,[',''' num2str(seed) '''']);
                string = AddLine(string,[',''' out_stopflag '''']);
                string = AddLine(string,[');']);
                res = exec(dbConn, string);
                
                %insert bestever_x
                for k = 1 : size(out_bestever_x,2)
                    string = ['INSERT INTO `' indbname '`.`bestever_x` (`id`, `setting_id`'];
                    string = AddLine(string,',`x` ');
                    string = AddLine(string,') VALUES (NULL');
                    string = AddLine(string,[',''' num2str(id) '''']);
                    string = AddLine(string,[',''' num2str(out_bestever_x(k)) '''']);
                    string = AddLine(string,[');']);
                    res = exec(dbConn, string);
                end
                
                %insert xstart
                for k = 1 : size(out_xstart,2)
                    string = ['INSERT INTO `' indbname '`.`xstart` (`id`, `setting_id`'];
                    string = AddLine(string,',`x` ');
                    string = AddLine(string,') VALUES (NULL');
                    string = AddLine(string,[',''' num2str(id) '''']);
                    string = AddLine(string,[',''' num2str(out_xstart(k)) '''']);
                    string = AddLine(string,[');']);
                    res = exec(dbConn, string);
                end
                
                %insert weights
%                 for k = 1 : size(out_weights,2)
%                     string = ['INSERT INTO `' indbname '`.`weights` (`id`, `setting_id`'];
%                     string = AddLine(string,',`weight` ');
%                     string = AddLine(string,') VALUES (NULL');
%                     string = AddLine(string,[',''' num2str(id) '''']);
%                     string = AddLine(string,[',''' num2str(out_weights(k)) '''']);
%                     string = AddLine(string,[');']);
%                     res = exec(dbConn, string);
%                 end
                
                cd ('..');
                [status, message, messageid] = rmdir(['output' num2str(i) 'repeat' num2str(j)],'s');
                disp(['output' num2str(i) 'repeat' num2str(j) ' added to DB']);
                end
            end
        end
        
        
        
        
        
        
        close(dbConn); 
        else
            disp(sprintf('Connection failed: %s', dbConn.Message));
        end

        % Close the connection so we don't run out of MySQL threads
       
     
end