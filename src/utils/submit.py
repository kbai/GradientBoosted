import mechanize
import sys
print sys.argv[1]
br = mechanize.Browser()
br.open("http://cs156.caltech.edu/scoreboard/")
br.select_form(nr=0)
br.form['teamid'] = 'jtqteaxb'
br.form.add_file(open(sys.argv[1]),'text/plain',sys.argv[1])
#br.form['file'] = './aaa.txt'
#br.form['valset'] = 1
req = br.submit()
print req.read()

