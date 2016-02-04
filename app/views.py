__author__ = 'liwang'

from app import app
import pandas as pd
import pymysql
#app = Flask(__name__)
from flask import render_template, request

from get_drug_list import get_drug_list
from plot_result import plot_result


@app.route('/')
@app.route('/index')
def index():
    mydb = pymysql.connect(host='localhost',user='root',password='wxx',db='Insight')

    query="SELECT * FROM Insight.CELL_LINE_INFO;"
    with mydb:
        cur = mydb.cursor()
        #just select the city from the world_innodb that the user inputs
        cur.execute(query)
        query_results = cur.fetchall()

    df_cell_info=pd.DataFrame(list(query_results),columns=['idx','cell_line_name','Primary_site','Tumour_origin','Comments'])
    cell_info_dict = get_drug_list(df_cell_info)

    return render_template("index.html",cell_info_dict=cell_info_dict)

@app.route('/about')
def pains_train_about():
    return render_template("about.html")

@app.route('/contact')
def pains_train_contact():
    return render_template("contact.html")

@app.route('/input')
def smile_input():
    return render_template("input.html")


@app.route('/output')
def synergize_output():
    cell_line_nm = request.args.get("lname")
    mydb = pymysql.connect(host='localhost',user='root',password='wxx',db='Insight')

    query="SELECT * FROM Insight.Prediction;"
    with mydb:
        cur = mydb.cursor()
        cur.execute(query)
        query_results = cur.fetchall()
    df_pred=pd.DataFrame(list(query_results),columns=['idx','cell_line_name','Combination_ID','Prob'])

    query="SELECT * FROM Insight.Drug_info;"

    with mydb:
        cur = mydb.cursor()
        cur.execute(query)
        query_results2 = cur.fetchall()
    df_drug=pd.DataFrame(list(query_results2),columns=['idx','DrugID','Target_genes'])
    df_match=df_pred.loc[df_pred['cell_line_name']==cell_line_nm,:]

    plotdiv=plot_result(df_match,df_drug,cell_line_nm)
    return render_template("output.html", plot=plotdiv)


@app.route('/error')
def error_output():
    return render_template("error.html")


if __name__ == "__main__":
    app.run(host='0.0.0.0', port=5000, debug=True)
