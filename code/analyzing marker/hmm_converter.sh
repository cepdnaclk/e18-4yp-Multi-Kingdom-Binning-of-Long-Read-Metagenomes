for i in *.hmm; do
  # Extract base filename without extension
  base_name="${i%.*}"

  # Perform HMM conversion and capture output in a variable
  converted_data=$(hmmconvert "$i")

  # Check if conversion was successful (exit code 0 indicates success)
  if [ $? -eq 0 ]; then
    # Write converted data to new file with base name and "_new.hmm" suffix
    echo "$converted_data" > "${base_name}_new.hmm"
    # Remove the original file after successful conversion
    $(rm "$i")
  else
    echo "Error: Conversion failed for file '$i'"
  fi
done
