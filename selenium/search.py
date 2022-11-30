from selenium import webdriver
from selenium.webdriver.common.by import By

driver = webdriver.Chrome()

driver.get("https://www.puzzle-minesweeper.com/mosaic-5x5-easy/")


elements = driver.find_elements(By.CLASS_NAME,"number")


for x in elements:
    if(x.text == ""):
        print(-1)
    else:
        print(x.text)
    
