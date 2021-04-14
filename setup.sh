mkdir -p ~/.streamlit/

echo "\
[general]\n\
email = \"aakanksha.gubbala@gmail.com\"\n\
" > ~/.streamlit/credentials.toml

echo "\
[theme]
primaryColor = '#0013F7'
backgroundColor = '#FFFFFF'
secondaryBackgroundColor = '#12E279'
textColor = '#000000'
font = 'monospace'
[server]\n\
headless = true\n\
enableCORS=false\n\
port = $PORT\n\
" > ~/.streamlit/config.toml