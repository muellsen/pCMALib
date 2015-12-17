function [dbConn] = DB_connect(open)

if (open)
% Database Server
host = 'localhost';

% Database Username/Password
user = 'root';
password = '';

% Database Name
dbName = 'CEC2005'; 

% JDBC Parameters
jdbcString = sprintf('jdbc:mysql://%s/%s', host, dbName);
jdbcDriver = 'com.mysql.jdbc.Driver';

% Set this to the path to your MySQL Connector/J JAR
javaaddpath('C:\Users\Georg\Downloads\mysql-connector-java-5.1.10\mysql-connector-java-5.1.10\mysql-connector-java-5.1.10-bin.jar')

% Create the database connection object
dbConn = database(dbName, user , password, jdbcDriver, jdbcString);


% Check to make sure that we successfully connected
if isconnection(dbConn)
% If the connection failed, print the error message
else
    disp(sprintf('Connection failed: %s', dbConn.Message));
end

% Close the connection so we don't run out of MySQL threads
close(dbConn); 
else
   
end
   

