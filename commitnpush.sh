cd 'C:\Users\Colin\Desktop\Master thesis\projet R\git\Phase-type'

git add .

echo 'Enter the commit message:'
read commitMessage

git commit -m "$commitMessage"

echo 'Enter the name of the branch:'
read branch

git push origin $branch

read