import requests
from bs4 import BeautifulSoup

# URL of the webpage containing the table
url = 'https://www.deciphergenomics.org/disorders/syndromes/list'

# Send a GET request to the webpage
response = requests.get(url)

# Check if the request was successful
if response.status_code == 200:
    # Parse the HTML content
    soup = BeautifulSoup(response.content, 'html.parser')

    # Find the table element
    table = soup.find('table')

    print(soup)
    # Extract table rows
    rows = table.find_all('tr')

    # Open a text file to write the scraped data
    with open('scraped_data.txt', 'w') as f:
        # Iterate over rows and write data to the file
        for row in rows:
            # Extract table cells (columns)
            cells = row.find_all(['th', 'td'])
            # Write data from each cell to the file
            for cell in cells:
                f.write(cell.get_text(strip=True) + '\t')
            f.write('\n')
    print('Scraping and exporting completed.')
else:
    print('Failed to retrieve webpage:', response.status_code)
